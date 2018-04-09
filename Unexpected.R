library(preprocessCore)
library(ggplot2)
library(dplyr)
library(reshape2)
library(cowplot)
library(Hmisc)

# Data identifier that prefixes filenames
data_name <- "Igen"

# Is Expression data already logged?
logged <- T

# Read and normalize data ----
# Set working directory and load data, omitting genes with missing data
setwd(paste("Z:/Research/", data_name, sep = ""))

Class_file <- read.table(paste(data_name, "_Class.txt", sep = ""), sep = "\t", header = T)
temp <- read.table(paste(data_name, "_Expression5.txt", sep = ""), sep = "\t", header = T)
Expression_file <- na.omit(temp)
Express <- as.matrix.data.frame(Expression_file)

# Define classes (states), sample sizes and number of genes
C0 <- Class_file$Class == 0
C1 <- Class_file$Class == 1
n0 <- sum(C0)
n1 <- sum(C1)
g <- dim(Express)[1]
c_g <- log(g) - digamma(1) + 1 / (2 * g)

# Quantile normalize the data for use with t-based method
if (logged) Express <- 2 ^ Express
Express[, C0] <- log2(normalize.quantiles(Express[, C0]))
Express[, C1] <- log2(normalize.quantiles(Express[, C1]))

# Calculate t and ordering statistics ----
diff <- apply(Express[, C0], 1, mean) - apply(Express[, C1], 1, mean)
row_vars <- cbind(apply(Express[, C0], 1, var), apply(Express[, C1], 1, var))
denom <- sqrt(row_vars[, 1] / n0 + row_vars[, 2] / n1)
t_stat <- diff / denom
nu <- denom ^ 4 / (row_vars[, 1] ^ 2 / (n0 ^ 2 * (n0 - 1)) + 
                      row_vars[, 2] ^ 2 / (n1 ^ 2 * (n1 - 1)))
p_val <- 2 * pt(t_stat, df = nu)
for (i in 1:g) if (p_val[i] > 1) p_val[i] <- 2 - p_val[i]

# Check for issues
subset(p_val, is.na(p_val)==TRUE)

# Add t-based statistics to observational data and order ----
working <- cbind(Express, p_val, diff)
row.names(working) <- row.names(Expression_file)
worked <- as.data.frame(working[order(p_val), ])
worked$p_rank <- 1:g
worked$adj_p_val <- g * c_g * worked$p_val / worked$p_rank

# Create data heat plot and Simes plot
Exp_img <- melt(Express)
Exp_img[, 1] <- kronecker(rep(1, n0 + n1), 1:g)
Data_plot <- ggplot(Exp_img, aes(Var1, Var2, fill=value)) + xlab("genes") + ylab("Subjects") + 
  geom_raster(show.legend = F) + geom_hline(yintercept = 20.45) +
  scale_fill_gradientn(colours = grey(seq(0.4,1,l=10))) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), 
        axis.text.y = element_blank()) +
  coord_cartesian(xlim = c(0, g), ylim = c(0, n1 + n0), expand = F)
Simes_plot <- ggplot(data = worked, aes(x = p_rank, y = p_val)) + xlab("gene rank") + 
  ylab(expression(p[(v)]^(1))) + geom_point(shape = ".") + 
  geom_path(aes(y = 0.05 * p_rank / (g * c_g)), linetype = "solid") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
plot_grid(Data_plot, Simes_plot, labels = "AUTO", vjust = 1, scale = c(1., 1.))

# Find set of differentially expressed genes for t-tests and write data to file ----
M <- worked %>% filter(adj_p_val < 0.05) %>% summarise(M = max(p_rank))
m <- M[1, 1]

write.table(worked[1:m, ], paste(data_name, "_diff_exped_t.txt", sep = ""))

# Calculate and classify correlations for t-based methods ----
# Declare variables
pairs <- choose(m, 2)
puc_dat <- matrix(NA, 50, 7)
FKG_vect <- matrix(NA, pairs, 1)
p0_vect <- matrix(NA, pairs, 1)
p1_vect <- matrix(NA, pairs, 1)
C0 <- c(C0, rep(FALSE, dim(worked)[2] - n0 - n1))
C1 <- c(C1, rep(FALSE, dim(worked)[2] - n0 - n1))

# Define FKG function
classify_FKG <- function(n0, n1, i, j) {
  cor_0 <- cor(twed[C0, i], twed[C0, j + 1])
  cor_1 <- cor(twed[C1, i], twed[C1, j + 1])
  diff_prod <- twed[(n1 + n0 + 2), i] * twed[(n1 + n0 + 2), j + 1]
  p_0 <- 2 * ifelse(cor_0 > 0, 1 - pt(cor_0 * sqrt((n0 - 2) / (1 - cor_0 ^ 2)), n0 - 2),
                    pt(cor_0 * sqrt((n0 - 2) / (1 - cor_0 ^ 2)), n0 - 2))
  p_1 <- 2 * ifelse(cor_1 > 0, 1 - pt(cor_1 * sqrt((n1 - 2) / (1 - cor_1 ^ 2)), n1 - 2),
                    pt(cor_1 * sqrt((n1 - 2) / (1 - cor_1 ^ 2)), n1 - 2))
  class_0 <- ifelse(cor_0 > 0, ifelse(diff_prod > 0, 1, 4), ifelse(diff_prod > 0, 2, 3))
  class_1 <- ifelse(cor_1 > 0, ifelse(diff_prod > 0, 1, 4), ifelse(diff_prod > 0, 2, 3))
  FKG_class <- ifelse(class_0 == class_1, class_0, 5)
  fn_out <- c(FKG_class, p_0, p_1)
  return(fn_out)
}

twed <- t(worked)  # transposing data for ease of use in function

# Use function to obtain correlations and p-values from t-tests and whether expected
ptm <- proc.time()
k <- 1
for (i in 1:(m - 1)) {
  for (j in i:(m - 1)) {
    temp <- classify_FKG(n0, n1, i, j)
    FKG_vect[k] <- temp[1]
    p0_vect[k] <- temp[2]
    p1_vect[k] <- temp[3]
    k <- k + 1
  }
  if(i %% 100 == 0) print(1 - (1 - i/m)*(1 - i/(m-1)))
}
proc.time() - ptm

# Write correlation data to file
write.table(cbind(FKG_vect, p0_vect, p1_vect), paste(data_name, "_corr_data_t.txt", sep = ""))
print("Correlation data table written")

# Create p-value comparison vectors
ptm <- proc.time()
c_pairs <- log(pairs) - digamma(1) + 1 / (2 * pairs)
MC_p0 <- pairs * c_pairs * p0_vect / rank(p0_vect)
print("MC state 0 p-values scaled for testing")
MC_p1 <- pairs * c_pairs * p1_vect / rank(p1_vect)
print("MC state 1 p-values scaled for testing")
proc.time() - ptm

# Summarize data for plotting
ptm <- proc.time()
for (n in 1:50) {
  sig_cor_tab <- NULL
  for (i in 1:5) {
    sig_cor_tab[i] <- length(which(subset(FKG_vect, MC_p0 < n * 0.02 & MC_p1 < n * 0.02) == i))
  }
  puc_dat[n, ] <- c(n * 0.02, 1 - (sig_cor_tab[1] + sig_cor_tab[3])/sum(sig_cor_tab), sig_cor_tab/sum(sig_cor_tab))
  print(paste("puc_dat is ", 2 * n, "% complete.", sep = ""))
}
proc.time() - ptm
write.table(puc_dat, paste(data_name, "_puc_data_t.txt", sep = ""))

# Define function used to create Proportion of Type figures
make.plot <- function(puc_dat){
  colnames(puc_dat) <- c("FDR", "AU", "PE", "PU", "NE", "NU", "MU")
  puc_dat.df <- melt(as.data.frame(puc_dat), id.vars = "FDR")
  temp <- ggplot(puc_dat.df, aes(FDR, value, group = variable)) + 
    geom_point(aes(shape = variable)) + theme_bw() +
    scale_shape_manual(values = c(3, 17, 2, 15, 0, 5)) +
    xlab(expression(paste(alpha, "'"))) + 
    ylab("Proportion of Type") + ylim(c(0, 1)) + 
    theme(legend.position="none", text = element_text(size = 16)) +
    coord_cartesian(xlim = c(-.02, 1.02), ylim = c(-.02, 1.02), expand = F)
  return(temp)
}

# Create summary plot for t-based methods
p1 <- make.plot(puc_dat)


# Rank normalize the data for use with normal approximation rank method
C0 <- Class_file$Class == 0
C1 <- Class_file$Class == 1
Express <- as.matrix.data.frame(Expression_file)
Express[, C0] <- apply(Express[, C0], 2, rank) / g
Express[, C1] <- apply(Express[, C1], 2, rank) / g

# Initialize variables
n_min <- min(n1, n0)
diff_n <- rep(NA, g)
p_val_n <- rep(NA, g)
C_min <- C0
if(n1 < n0) C_min <- C1
mu_n <- n1 * n0 / 2
sigma_n <- sqrt(n1 * n0 * (n1 + n0 + 1) / 12)

# Calculate differences and p-vals for normal approx
diff_n <- apply(Express[, C0], 1, mean) - apply(Express[, C1], 1, mean)
rank_sum <- apply(apply(Express, 1, rank)[C_min, ], 2, sum)
U_min <- rank_sum - n_min * (n_min + 1) / 2
Z_n <- abs((U_min - mu_n) / sigma_n)
p_val_n <- 2 * pnorm(Z_n, lower.tail = F)

# Check for issues
subset(p_val_n, is.na(p_val_n)==TRUE)

# Add rank based statistics to observational data and order ----
working_n <- cbind(Express, p_val_n, diff_n)
row.names(working_n) <- row.names(Expression_file)
worked_n <- as.data.frame(working_n[order(p_val_n), ])
worked_n$p_rank <- 1:g
worked_n$adj_p_val <- g * c_g * worked_n$p_val_n / worked_n$p_rank

# Create Simes plot
ggplot(data = worked_n, aes(x = p_rank, y = p_val_n)) + xlab("gene rank") + 
  ylab("adjusted p-values") + geom_point(shape = ".") + 
  geom_path(aes(y = 0.05 * p_rank / (g * c_g)), col = "red") +
  theme_bw() + coord_cartesian(xlim = c(0, g), ylim = c(0, 1), expand = F)

# Find set of differentially expressed genes for t-tests and write data to file ----
M_n <- worked_n %>% filter(adj_p_val < 0.05) %>% summarise(M = max(p_rank))
m_n <- M_n[1, 1]

write.table(worked_n[1:m_n, ], paste(data_name, "_diff_exped_n.txt", sep = ""))

# Calculate and classify correlations for t-based methods ----
# Declare variables
pairs <- choose(m_n, 2)
puc_dat <- matrix(NA, 50, 7)
FKG_vect <- matrix(NA, pairs, 1)
p0_vect <- matrix(NA, pairs, 1)
p1_vect <- matrix(NA, pairs, 1)
C0 <- c(C0, rep(FALSE, dim(worked_n)[2] - n0 - n1))
C1 <- c(C1, rep(FALSE, dim(worked_n)[2] - n0 - n1))

# Define FKG function for rank correlations
classify_FKG_n <- function(n0, n1, i, j) {
  temp1 <- rcorr(twed[C0, i], twed[C0, j + 1], type = "spearman")
  temp2 <- rcorr(twed[C1, i], twed[C1, j + 1], type = "spearman")
  cor_0 <- temp1$r[1, 2]
  cor_1 <- temp2$r[1, 2]
  diff_prod <- worked_n$diff[i] * worked_n$diff[j + 1]
  class_0 <- ifelse(cor_0 > 0, ifelse(diff_prod > 0, 1, 4), ifelse(diff_prod > 0, 2, 3))
  class_1 <- ifelse(cor_1 > 0, ifelse(diff_prod > 0, 1, 4), ifelse(diff_prod > 0, 2, 3))
  FKG_class <- ifelse(class_0 == class_1, class_0, 5)
  fn_out <- c(FKG_class, temp1$P[1, 2], temp2$P[1, 2])
  return(fn_out)
}

twed <- t(worked_n)  # transposing data for ease of use in function

# Use function to obtain correlations and p-values from rank tests and whether expected
ptm <- proc.time()
k <- 1
for (i in 1:(m - 1)) {
  for (j in i:(m - 1)) {
    temp <- classify_FKG_n(n0, n1, i, j)
    FKG_vect[k] <- as.integer(temp[1])
    p0_vect[k] <- temp[2]
    p1_vect[k] <- temp[3]
    k <- k + 1
  }
  if(i %% 100 == 0) print(1 - (1 - i/m)*(1 - i/(m-1)))
}
proc.time() - ptm


# Write correlation data to file
ptm <- proc.time()
write.table(cbind(FKG_vect, p0_vect, p1_vect), paste(data_name, "_corr_data_n.txt", sep = ""))
print("Correlation data table written")
proc.time() - ptm

# Create p-value comparison vectors
ptm <- proc.time()
c_pairs <- log(pairs) - digamma(1) + 1 / (2 * pairs)
MC_p0 <- pairs * c_pairs * p0_vect / rank(p0_vect)
print("MC state 0 p-values scaled for testing")
MC_p1 <- pairs * c_pairs * p1_vect / rank(p1_vect)
print("MC state 1 p-values scaled for testing")
proc.time() - ptm

# Summarize data for plotting
ptm <- proc.time()
for (n in 1:50) {
  sig_cor_tab <- NULL
  for (i in 1:5) {
    sig_cor_tab[i] <- length(which(subset(FKG_vect, MC_p0 < n * 0.02 & MC_p1 < n * 0.02) == i))
  }
  puc_dat[n, ] <- c(n * 0.02, 1 - (sig_cor_tab[1] + sig_cor_tab[3])/sum(sig_cor_tab), sig_cor_tab/sum(sig_cor_tab))
  print(paste("puc_dat is ", 2 * n, "% complete.", sep = ""))
}
proc.time() - ptm
write.table(puc_dat, paste(data_name, "_puc_data_n.txt", sep = ""))

# Create summary plot for rank based methods
p2 <- make.plot(puc_dat)

# Create summary figure
plot_grid(p1, p2, labels = "AUTO", hjust = 0, vjust = 1, scale = c(1., 1.))


# Read in specific rank distribution
pval_tab <- read.table(paste(data_name, "_pval_tab.txt", sep = ""))

# Initialize variables
n_min <- min(n1, n0)
med_sum <- (n1 + n0 + 1) * n_min / 2
len <- sum(as.numeric(pval_tab$Freq))
min_sum <- sum(1:n_min)
len_tab <- length(pval_tab$Freq)

diff_r <- rep(NA, g)
rank_sum <- rep(NA, g)
p_val_r <- rep(NA, g)
pre_p_val <- rep(NA, g)
rand_p_val <- rep(NA, g)
C0 <- Class_file$Class == 0
C1 <- Class_file$Class == 1
C_min <- C0
if(n1 < n0) C_min <- C1

# find p-vals with use of null distribution
for (i in 1:g) {
  diff_r[i] <- mean(Express[i, C0]) - mean(Express[i, C1])
  rank_sum[i] <- sum(rank(Express[i, ])[C_min])
  opp_sum <- 2 * med_sum - rank_sum[i]
  a <- sum(as.numeric(pval_tab$Freq[rank_sum[i] < pval_tab$Stat]))
  b <- sum(as.numeric(pval_tab$Freq[rank_sum[i] <= pval_tab$Stat]))
  c <- sum(as.numeric(pval_tab$Freq[opp_sum < pval_tab$Stat]))
  d <- sum(as.numeric(pval_tab$Freq[opp_sum <= pval_tab$Stat]))
  p_0 <- 1 + (min(a, d) - max(a, d)) / len
  p_1 <- 1 + (min(b, c) - max(b, c)) / len
  pre_p_val[i] <- min(p_0, p_1)
  p_val_r[i] <- max(p_0, p_1)
  if(p_val_r[i] > 1) p_val_r[i] <- 1
}

subset(p_val_r, is.na(p_val_r)==TRUE)

# sim randomized p-values and find rejection prob
rej_sum <- rep(0, g)
nsim <- 10000
p_store <- matrix(NA, g, nsim)
c_g <- log(g) - digamma(1) + 0.5 / g
for (i in 1:nsim) {
  rand_p_val <- runif(g, pre_p_val, p_val_r)
  p_store[, i] <- rand_p_val
  p_rank <- rank(rand_p_val)
  test_val <- g * c_g * rand_p_val / p_rank
  df <- data.frame(p_rank, test_val)
  M <- df %>% filter(test_val < 0.05) %>% summarise(M = max(p_rank))
  rej_sum <- rej_sum + as.integer(p_rank <= M[1, 1])
}
rej_prob <- rej_sum / nsim

rej_prob <- rej_prob[order(p_val_r)]
p_store <- p_store[order(p_val_r), ]

m_r <- sum(rej_prob >= 0.85)

# Add simulated statistics to observational data and order ----
working_r <- cbind(Express, p_val_r, diff_r)
row.names(working_r) <- row.names(Expression_file)
worked_r <- as.data.frame(working_r[order(p_val_r), ])
worked_r$p_rank <- 1:g
worked_r$adj_p_val <- g * c_g * worked_r$p_val_r / worked_r$p_rank

# Find set of differentially expressed genes for rank tests and write data to file ----
M_s <- worked_r %>% filter(adj_p_val < 0.05) %>% summarise(M = max(p_rank))
m_s <- M[1, 1]
worked_s <- worked_r[1:m_s, ]
worked_r <- worked_r[1:m_r, ]

write.table(worked_s, paste(data_name, "_diff_exped_s.txt", sep = ""))
write.table(worked_r, paste(data_name, "_diff_exped_r.txt", sep = ""))

# Read in specific simulated correlation distributions
n1_null <- read.table(paste(data_name, "_n1_null_fin.txt", sep = ""))
n2_null <- read.table(paste(data_name, "_n2_null_fin.txt", sep = ""))
cnt_scale <- sum(n1_null$Freq)

# Initialize correlation variables
pairs_s <- choose(m_s, 2)
pairs_r <- choose(m_r, 2)
puc_dat_s <- matrix(NA, 50, 7)
puc_dat_r <- matrix(NA, 50, 7)
rii_s <- matrix(NA, pairs_s, 1)
rii_r <- matrix(NA, pairs_r, 1)
p0_vect_s <- matrix(NA, pairs_s, 1)
p1_vect_s <- matrix(NA, pairs_s, 1)
p0_vect_r <- matrix(NA, pairs_r, 1)
p1_vect_r <- matrix(NA, pairs_r, 1)

# Define randomized FKG function
classify_FKG_r <- function(i, j) {
  temp1 <- rcorr(twed_r[C0, i], twed_r[C0, j + 1], type = "spearman")
  temp2 <- rcorr(twed_r[C1, i], twed_r[C1, j + 1], type = "spearman")
  cor_0 <- temp1$r[1, 2]
  cor_1 <- temp2$r[1, 2]
  diff_prod <- worked_r$diff[i] * worked_r$diff[j + 1]
  p_0 <- 0.5 * sum(n1_null$Freq[n1_null$Stat >= abs(cor_0)]) / cnt_scale
  pre_p_0 <- 0.5 * sum(n1_null$Freq[n1_null$Stat > abs(cor_0)]) / cnt_scale
  ran_p_0 <- ifelse(p_0 == pre_p_0, p_0, runif(1, pre_p_0, p_0))
  p_1 <- 0.5 * sum(n2_null$Freq[n2_null$Stat >= abs(cor_1)]) / cnt_scale
  pre_p_1 <- 0.5 * sum(n2_null$Freq[n2_null$Stat > abs(cor_1)]) / cnt_scale
  ran_p_1 <- ifelse(p_1 == pre_p_1, p_1, runif(1, pre_p_1, p_1))
  class_0 <- ifelse(cor_0 > 0, ifelse(diff_prod > 0, 1, 4), ifelse(diff_prod > 0, 2, 3))
  class_1 <- ifelse(cor_1 > 0, ifelse(diff_prod > 0, 1, 4), ifelse(diff_prod > 0, 2, 3))
  FKG_class <- ifelse(class_0 == class_1, class_0, 5)
  fn_out <- c(FKG_class, ran_p_0, ran_p_1)
  return(fn_out)
}

# Define simulated FKG function
classify_FKG_s <- function(i, j) {
  temp1 <- rcorr(twed_s[C0, i], twed_s[C0, j + 1], type = "spearman")
  temp2 <- rcorr(twed_s[C1, i], twed_s[C1, j + 1], type = "spearman")
  cor_0 <- temp1$r[1, 2]
  cor_1 <- temp2$r[1, 2]
  diff_prod <- worked_s$diff[i] * worked_s$diff[j + 1]
  p_0 <- 0.5 * sum(n1_null$Freq[n1_null$Stat >= abs(cor_0)]) / cnt_scale
  p_1 <- 0.5 * sum(n2_null$Freq[n2_null$Stat >= abs(cor_1)]) / cnt_scale
  class_0 <- ifelse(cor_0 > 0, ifelse(diff_prod > 0, 1, 4), ifelse(diff_prod > 0, 2, 3))
  class_1 <- ifelse(cor_1 > 0, ifelse(diff_prod > 0, 1, 4), ifelse(diff_prod > 0, 2, 3))
  FKG_class <- ifelse(class_0 == class_1, class_0, 5)
  fn_out <- c(FKG_class, p_0, p_1)
  return(fn_out)
}

# transposing data for ease of use in function
twed_s <- t(worked_s)
twed_r <- t(worked_r)

# Use function to obtain correlations and p-values from simulated rank tests and whether expected
ptm <- proc.time()
k <- 1
for (i in 1:(m_s - 1)) {
  for (j in i:(m_s - 1)) {
    temp <- classify_FKG_s(i, j)
    rii_s[k] <- as.integer(temp[1])
    p0_vect_s[k] <- temp[2]
    p1_vect_s[k] <- temp[3]
    k <- k + 1
  }
  if(i %% 100 == 0) print(1 - (1 - i/m_s) * (1 - i/(m_s - 1)))
}
proc.time() - ptm

# Write simulated correlation data to file
ptm <- proc.time()
write.table(cbind(rii_s, p0_vect_s, p1_vect_s), paste(data_name, "_corr_data_s.txt", sep = ""))
print("simulated corr table written")
proc.time() - ptm

# Use function to obtain correlations and p-values from randomized rank tests and whether expected
ptm <- proc.time()
k <- 1
for (i in 1:(m_r - 1)) {
  for (j in i:(m_r - 1)) {
    temp <- classify_FKG_r(i, j)
    rii_r[k] <- as.integer(temp[1])
    p0_vect_r[k] <- temp[2]
    p1_vect_r[k] <- temp[3]
    k <- k + 1
  }
  if(i %% 100 == 0) print(1 - (1 - i/m_r) * (1 - i/(m_r - 1)))
}
proc.time() - ptm

# Write randomized correlation data to file
ptm <- proc.time()
write.table(cbind(rii_r, p0_vect_r, p1_vect_r), paste(data_name, "_corr_data_r.txt", sep = ""))
print("randomized corr table written")
proc.time() - ptm

# Create p-value comparison vectors
ptm <- proc.time()
c_pairs_s <- log(pairs_s) - digamma(1) + 1 / (2 * pairs_s)
MC_p0_s <- pairs_s * c_pairs_s * p0_vect_s / rank(p0_vect_s)
print("MC state 0 simulated p-values scaled for testing")
MC_p1_s <- pairs_s * c_pairs_s * p1_vect_s / rank(p1_vect_s)
print("MC state 1 simulated p-values scaled for testing")

c_pairs_r <- log(pairs_r) - digamma(1) + 1 / (2 * pairs_r)
MC_p0_r <- pairs_r * c_pairs_r * p0_vect_r / rank(p0_vect_r)
print("MC state 0 randomized p-values scaled for testing")
MC_p1_r <- pairs_r * c_pairs_r * p1_vect_r / rank(p1_vect_r)
print("MC state 1 randomized p-values scaled for testing")
proc.time() - ptm

# Summarize data for plotting
ptm <- proc.time()
for (n in 1:50) {
  sig_cor_tab <- NULL
  for (i in 1:5) {
    sig_cor_tab[i] <- length(which(subset(rii_r, MC_p0_r < n * 0.02 & MC_p1_r < n * 0.02) == i))
  }
  puc_dat_r[n, ] <- c(n * 0.02, 1 - (sig_cor_tab[1] + sig_cor_tab[3])/sum(sig_cor_tab), sig_cor_tab/sum(sig_cor_tab))
  print(paste("puc_dat_r is ", 2 * n, "% complete.", sep = ""))
}
proc.time() - ptm
write.table(puc_dat_r, paste(data_name, "_puc_data_r.txt", sep = ""))

ptm <- proc.time()
for (n in 1:100) {
  sig_cor_tab <- NULL
  for (i in 1:5) {
    sig_cor_tab[i] <- length(which(subset(rii_s, MC_p0_s < n * 0.02 & MC_p1_s < n * 0.02) == i))
  }
  puc_dat_s[n, ] <- c(n * 0.02, 1 - (sig_cor_tab[1] + sig_cor_tab[3])/sum(sig_cor_tab), sig_cor_tab/sum(sig_cor_tab))
  print(paste("puc_dat_s is ", 2 * n, "% complete.", sep = ""))
}
proc.time() - ptm
write.table(puc_dat_s, paste(data_name, "_puc_data_s.txt", sep = ""))

# Create summary plot for rank based methods
p3 <- make.plot(puc_dat_s)
p4 <- make.plot(puc_dat_r)

# Create summary figure
plot_grid(p3, p4, labels = "AUTO", hjust = 0, vjust = 1, scale = c(1., 1.))