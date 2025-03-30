rm(list=ls(all=TRUE))
library(MCMCpack)
library(mvtnorm)
library(rstudioapi)
library(ggplot2)
library(scales)

script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
script_path <- getSourceEditorContext()$path
script_dir <- dirname(script_path)

file_path <- file.path(script_dir, "GoogleSearchIndex.txt")
dat <- read.delim(file_path)

dat$Week <- as.Date(dat$Week)

# EDA
ggplot(dat, aes(x = Week, y = covid)) +
  geom_line(color = "brown1", linewidth = 1) +  
  geom_point(color = "brown1", size = 3) +    
  labs(
    title = "Public Interest in COVID-19 Over Time (Google Search Index)",
    x = "Date",
    y = "Search Index"
  ) +
  theme_minimal() +  
  scale_x_date(
    date_breaks = "3 months",       
    date_labels = "%Y-%m",          
    limits = c(min(dat$Week), max(dat$Week))     
  )


y=as.numeric(dat$covid)
print(length(y))


#### Model setup
n.all = n.all = length(y)
p=1 ## order of AR process
K=2 ## number of mixing component

Y = matrix(y[(p + 1):n.all], ncol = 1) ## y_{p+1:T}

Fmtx = sapply(1:p, function(k) {
  y[(p + 1 - k):(n.all - k)]
})
Fmtx = t(Fmtx)
print(Fmtx)

n=length(Y) ## T-p 

## prior hyperparameters
m0 = matrix(rep(0, p), ncol=1)     ## AR 系数的先验均值（非信息性，0 向量）
C0 = 10 * diag(p)                  ## AR 系数先验协方差矩阵，较大表示不确定性高
n0 = 0.02                          ## 残差方差的 inverse-gamma 先验自由度
d0 = 0.02                          ## 残差方差的 inverse-gamma 先验尺度参数
a = rep(1, K)                      ## Dirichlet 先验参数（均匀先验），用于 omega 的采样


#使用的是 Dirichlet 共轭采样，后验参数为 a + 成分频数
sample_omega = function(L.cur){
  n.vec = sapply(1:K, function(k){ sum(L.cur == k) })  ## 统计每个成分被分配了多少数据点
  rdirichlet(1, a + n.vec)                            ## Dirichlet 后验采样权重
}


#对于单个时间点t，我们不知道它来自哪个成分k，这个函数尝试估计那个概率并采样
sample_L_one = function(beta.cur, omega.cur, nu.cur, y.cur, Fmtx.cur){
  prob_k = function(k){
    beta.use = beta.cur[((k-1)*p + 1):(k*p)]            ## 取出第 k 个 AR 成分的系数
    mu = sum(beta.use * Fmtx.cur)                       ## 当前时间点的预测值
    omega.cur[k] * (dnorm(y.cur, mean=mu, sd=sqrt(nu.cur)) + 1e-10)  ## 成分权重 × 对应的高斯密度值
  }

  prob.vec = sapply(1:K, prob_k)                         ## 得到每个成分的（非标准化）概率
  prob.vec = prob.vec + 1e-10

  L.sample = sample(1:K, 1, prob = (prob.vec + 1e-8) / sum(prob.vec + 1e-8)) ## 按照归一化后的概率进行采样
  return(L.sample)
}

#DS修改版本
sample_L_one <- function(beta.cur, omega.cur, nu.cur, y.cur, Fmtx.cur, DEBUG = TRUE) {
  prob_k <- function(k) {
    beta.use <- beta.cur[((k-1)*p + 1):(k*p)]            # 取出第 k 个 AR 成分的系数
    mu <- sum(beta.use * Fmtx.cur)                       # 当前时间点的预测值
    likelihood <- dnorm(y.cur, mean = mu, sd = sqrt(nu.cur))
    prob <- omega.cur[k] * (likelihood + 1e-10)          # 成分权重 × 高斯密度值
    
    if (DEBUG) {
      cat(sprintf(
        "[DEBUG prob_k] k=%d | beta.use=%s | mu=%.3f | likelihood=%.3e | omega=%.3f | prob=%.3e\n",
        k, paste(round(beta.use, 3), collapse = ","), mu, likelihood, omega.cur[k], prob
      ))
    }
    
    return(prob)
  }
  
  prob.vec <- sapply(1:K, prob_k)                       # 得到每个成分的非标准化概率
  prob.vec <- prob.vec + 1e-10                         # 避免数值下溢
  
  if (DEBUG) {
    cat("[DEBUG prob.vec] Raw prob.vec =", paste(round(prob.vec, 5), collapse = ", "), "\n")
    cat("[DEBUG prob.vec] Normalized prob.vec =", 
        paste(round(prob.vec / sum(prob.vec), 5), collapse = ", "), "\n")
  }
  
  L.sample <- sample(1:K, 1, prob = (prob.vec + 1e-8) / sum(prob.vec + 1e-8))  # 归一化后采样
  
  if (DEBUG) {
    cat("[DEBUG L.sample] Sampled L =", L.sample, "\n\n")
  }
  
  return(L.sample)
}


#采样整条序列的标签向量 𝐿1:𝑛
sample_L=function(y,x,beta.cur,omega.cur,nu.cur){
  L.new=sapply(1:n, function(j){sample_L_one(beta.cur,omega.cur,nu.cur,y.cur=y[j,],Fmtx.cur=x[,j])})
  return(L.new)
}

#采样噪声方差
sample_nu = function(L.cur, beta.cur) {
  n.star = n0 + n  #n.star为残差方差的 Inverse-Gamma 后验自由度参数
  err.y = function(idx) {
    L.use = L.cur[idx]                             # 当前观测属于哪个成分 k
    beta.use = beta.cur[((L.use - 1) * p + 1):(L.use * p)]  # 拿到该成分的 beta
    err = Y[idx, ] - sum(Fmtx[, idx] * beta.use)   # 计算该观测点的预测误差
    return(err^2)
  }
  d.star = d0 + sum(sapply(1:n, err.y))            # 后验 IG 的尺度参数
  1 / rgamma(1, shape = n.star / 2, rate = d.star / 2)  # 从 Inverse-Gamma 中采样
}

#采样成分k的AR系数
# 修改sample_beta函数中的部分代码
sample_beta = function(k, L.cur, nu.cur) {
  idx.select = (L.cur == k)
  n.k = sum(idx.select)
  if (n.k == 0) {
    m.k = m0
    C.k = C0
  } else {
    y.tilde.k = Y[idx.select, , drop = FALSE]  # 保持矩阵结构
    Fmtx.tilde.k = Fmtx[, idx.select, drop = FALSE]  # 保持矩阵结构
    e.k = y.tilde.k - t(Fmtx.tilde.k) %*% m0
    Q.k = t(Fmtx.tilde.k) %*% C0 %*% Fmtx.tilde.k + diag(n.k)
    Q.k.inv = chol2inv(chol(Q.k))
    A.k = C0 %*% Fmtx.tilde.k %*% Q.k.inv
    m.k = m0 + A.k %*% e.k
    C.k = C0 - A.k %*% Q.k %*% t(A.k)
  }
  
  rmvnorm(1, m.k, nu.cur * C.k)
}

# 确保在sample_L函数中传递Fmtx列时保持矩阵结构
sample_L = function(y, x, beta.cur, omega.cur, nu.cur) {
  L.new = sapply(1:n, function(j) {
    sample_L_one(beta.cur, omega.cur, nu.cur, y.cur = y[j, , drop = FALSE], Fmtx.cur = x[, j, drop = FALSE])
  })
  return(L.new)
}

# 总结，没有被分配到点的成分也必须更新，但是用的是先验


#### MCMC setup
nsim = 1000  ## 设定采样迭代次数

## 存储每次迭代后的参数结果
beta.mtx  = matrix(0, nrow = p*K, ncol = nsim)  ## 每行是 beta_1, ..., beta_K 的拼接
L.mtx     = matrix(0, nrow = n,    ncol = nsim) ## 每一列存储当前 L（哪个观测属于哪个成分）
omega.mtx = matrix(0, nrow = K,    ncol = nsim) ## 每次迭代混合权重向量 omega
nu.vec    = rep(0, nsim)                        ## 每次迭代噪声方差 nu

beta.cur  = rep(0, p*K)       ## 每个成分的 beta 系数初始设为 0
L.cur     = rep(1, n)         ## 初始假设所有点都属于第一个成分
omega.cur = rep(1/K, K)       ## 初始 omega 平均分配（均匀）
nu.cur    = 1                 ## 初始噪声方差设为 1


## Gibbs Sampler
# 设置测试模式开关
# 设置测试模式开关
TEST <- FALSE  # 设置为TRUE开启详细输出，FALSE关闭

for (i in 1:nsim) {
  set.seed(i)
  
  ## sample omega
  omega.cur <- sample_omega(L.cur)
  omega.mtx[, i] <- omega.cur
  
  # 测试模式：打印omega（修正括号）
  if (TEST) {
    print(paste("Iteration", i, "| omega =", 
                paste(round(omega.cur, 3), collapse=", ")))  # 注意此处闭合三个括号
  }
  
  ## sample L
  L.cur <- sample_L(Y, Fmtx, beta.cur, omega.cur, nu.cur)
  L.mtx[, i] <- L.cur
  
  # 测试模式：打印L的分布统计
  if (TEST) {
    l_counts <- table(factor(L.cur, levels=1:K))
    print(paste("Iteration", i, "| L counts:", paste(l_counts, collapse=", ")))
  }
  
  ## sample nu
  nu.cur <- sample_nu(L.cur, beta.cur)
  nu.vec[i] <- nu.cur
  
  # 测试模式：打印nu
  if (TEST) {
    print(paste("Iteration", i, "| nu =", round(nu.cur, 3)))
  }
  
  ## sample beta
  beta.cur <- as.vector(sapply(1:K, function(k) { sample_beta(k, L.cur, nu.cur) }))
  beta.mtx[, i] <- beta.cur
  
  # 测试模式：打印beta
  if (TEST) {
    print(paste("Iteration", i, "| beta =", 
                paste(round(beta.cur, 3), collapse=", ")))  # 同样修正此处括号
  }
  
  ## 常规进度提示（无论TEST是否开启）
  if (i %% 1000 == 0) {
    print(paste("Number of iterations:", i))
  }
}

#### show the result
#后 10,000 次被认为已经收敛到后验分布，可用于统计推断（如预测）
sample.select.idx=seq(5001,10000,by=1)

#根据后验样本进行预测
post.pred.y.mix=function(idx){
  
  k.vec.use=L.mtx[,idx]
  beta.use=beta.mtx[,idx]
  nu.use=nu.vec[idx]
  
  
  get.mean=function(s){
    k.use=k.vec.use[s]
    sum(Fmtx[,s]*beta.use[((k.use-1)*p+1):(k.use*p)])
  }
  mu.y=sapply(1:n, get.mean)
  ## 给第 k 个时间点，生成一个预测值
  sapply(1:length(mu.y), function(k){rnorm(1,mu.y[k],sqrt(nu.use))})
  
}  
#y.post.pred.sample是对所有点的后验预测
y.post.pred.sample=sapply(sample.select.idx, post.pred.y.mix)

summary.vec95=function(vec){
  c(unname(quantile(vec,0.025)),mean(vec),unname(quantile(vec,0.975)))
}

summary.y=apply(y.post.pred.sample,MARGIN=1,summary.vec95)

#plot posterior inference
library(ggplot2)

n <- length(Y)  # 有效观测长度
df_plot <- data.frame(
  Time   = 1:n,
  Truth  = Y,
  Mean   = summary.y[2, ],
  Lower  = summary.y[1, ],
  Upper  = summary.y[3, ]
)

# 下面用 ggplot2 绘图
p <- ggplot(df_plot, aes(x = Time)) +
  # 真实值：黑色线+点
  geom_line(aes(y = Truth, color = "Truth"), size = 0.8) +
  geom_point(aes(y = Truth, color = "Truth"), size = 1.5) +
  
  geom_line(aes(y = Mean, color = "Mean"), size = 0.8, linetype = "dashed") +
  geom_point(aes(y = Mean, color = "Mean"), size = 1.5, shape = 4) +
 
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = "95% C.I."), alpha = 0.2) +
  
  # color
  scale_color_manual(values = c("Truth" = "black", "Mean" = "grey")) +
  scale_fill_manual(values  = c("95% C.I." = "purple")) +
  
  labs(x = "Time", y = "", color = "", fill = "") +
  theme_minimal()

print(p)

check_normal_density <- function(x) {
  density <- dnorm(x, mean = 0, sd = 1)
  
  cat("Value:", x, "\n")
  cat("Density under N(0,1):", format(density, scientific = TRUE), "\n")
  
  if (density < 1e-10) {
    cat("⚠️  WARNING: This density is extremely small — potential underflow risk!\n")
  } else {
    cat("✅  Density is within normal range.\n")
  }
  
  return(invisible(density))
}
check_normal_density(81)      
check_normal_density(8) 


