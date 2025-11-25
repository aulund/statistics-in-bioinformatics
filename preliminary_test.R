birds <- c(3.2,3.7,4.2,5.2,2.8,4.3,2.9,3.8,2.2)

mean(birds)
median(birds)
shapiro.test(birds)

hist(birds)


snake <- c(182, 172, 121, 250, 224, 209)

confint.default(snake)
mean(snake)
median(snake)
shapiro.test(snake)

hist(snake)


t.test(snake)$conf.int


blood <- 1:9
bp <- c(140,137,135,131,129,125,124,120,118)
bmi <- c(40,41,35,31,28,29,22,25,21)
dat <- data.frame(blood, blood_pressure = bp, BMI = bmi)

cor(dat$blood_pressure, dat$BMI, method ="pearson")




pat <- 1:9
cig <- c(70, 75, 80, 82, 90, 91, 100, 105, 120)
age <- c(90,92,83,84,81,80,76,75,70)
cigdat <- data.frame(pat, Cigg = cig, Age = age)
fit <- lm(cigdat$Age ~ cigdat$Cigg)
summary(fit)
test <- 119.83 - 0.427 * 71

prostpat <- 1:7
old <- c(1.2, 2.4, 2.5, 4.0, 3.9, 0.9, 1.1)
young <- c(0.9, 0.5, 0.3, 0.9, 1.2, 1.3, 1.2)
prostdat <- data.frame(prostpat, Old = old, Young = young)

var.test(prostdat$Old, prostdat$Young)
t.test(prostdat$Old, prostdat$Young)

lake <- 1:3
a <- c(1.6, 0.2, 1.3, 2.9, 2.2)
b <- c(0.9, 0.8, 2.1, 0.3, 0.9)
c <- c(2.4, 1.2, 2.1, 0.3, 0.9)
y <- c(a, b , c)
lakes <- data.frame ( A=a, B=b, C=c)

fit <- aov(lakes$A,lakes$B,lakes$C)

lakes <- data.frame(A = a, B = b, C = c)

long <- stack(lakes)           # produces columns: values, ind
names(long) <- c("value", "group")

fit <- aov(value ~ group, data = long)
summary(fit)




# data
newborns <- c(120, 131, 100, 107, 108, 130, 99)
children <- c(80, 88, 81, 79, 129, 121, 115)

# quick summaries
n1 <- length(newborns); n2 <- length(children)
m1 <- mean(newborns); sd1 <- sd(newborns)
m2 <- mean(children); sd2 <- sd(children)
cat("Newborns: n=", n1, " mean=", round(m1,3), " sd=", round(sd1,3), "\n")
cat("Children: n=", n2, " mean=", round(m2,3), " sd=", round(sd2,3), "\n")

# F test for equality of variances
var.test(newborns, children)

# Two-sample t-tests
t.test(newborns, children, var.equal = TRUE)   # pooled t
t.test(newborns, children)                     # Welch (default)

# Nonparametric alternative
wilcox.test(newborns, children)

# Cohen's d (pooled)
sp <- sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2) / (n1 + n2 - 2))
d <- (m2 - m1) / sp    # Children - Newborns
d

# Create data frame for plotting
library(dplyr)
library(ggplot2)

dat <- data.frame(
  value = c(newborns, children),
  group = factor(rep(c("Newborns","Children"), each = n1))
)

summary_df <- dat %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mean = mean(value),
    sd = sd(value),
    se = sd / sqrt(n),
    ci = qt(0.975, df = n - 1) * se,
    ymin = mean - ci,
    ymax = mean + ci
  )

# bar chart with 95% CI error bars
ggplot(summary_df, aes(x = group, y = mean, fill = group)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.15, size = 0.7) +
  labs(x = "", y = "Mean value (with 95% CI)") +
  theme_minimal()

