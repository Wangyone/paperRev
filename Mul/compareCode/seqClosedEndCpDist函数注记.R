library(npcp)

# 1. 准备数据
# 假设 yt 是总数据，前 100 个是训练，后 50 个是监测
m <- 100 
train_data <- yt[1:m, , drop=FALSE] # 必须是矩阵
new_data   <- yt[(m+1):nrow(yt), , drop=FALSE] # 监测部分

# 关键：计算 n (总视窗)
# 监测期是 m+1 到 n
n <- m + nrow(new_data) 

# --- Step 1: 模拟 Null 分布 ---
# 注意：对于时间序列，必须用 method="mult"
# n 必须是 m + 监测长度
traj <- simClosedEndCpDist(x.learn = train_data, 
                           n = n,  # <--- 这里填总长度 150
                           method = "mult", 
                           B = 1000)

# --- Step 2: 计算阈值 ---
# alpha 控制误报率
thresh <- threshClosedEndCpDist(traj, alpha = 0.05)

# --- Step 3: 计算检测统计量 ---
# 这里传入 x.learn 和 新数据 x
# 这一步计算的是从 m+1 到 n 的所有统计量值
det <- detClosedEndCpDist(x.learn = train_data, 
                          x = new_data)

# --- Step 4: 实施监控 ---
# 比较 det 和 thresh
# 选择统计量类型，例如 "mac" (Tmq) 或 "mc" (Qmk)
res_mac <- monClosedEndCpDist(det, thresh, statistic = "mac")
res_mc  <- monClosedEndCpDist(det, thresh, statistic = "mc")

# --- 查看结果 ---
print(res_mac$alarm)      # 是否报警 (TRUE/FALSE)
print(res_mac$time.alarm) # 报警时刻 (相对于 m 的偏移量? 需检查文档)

# 文档细读：
# time.alarm: an integer corresponding to the time at which the detector function has exceeded...
# 注意：根据经验，这里的 time.alarm 通常是相对于 monitoring start 的索引 (1, 2...)
# 或者是绝对时间 (m+1, m+2...)，具体取决于包版本，但计算 EDD 时通常直接用这个值即可。