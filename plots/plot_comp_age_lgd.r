library(ggpubr)
library(gridExtra)
cand <- read.delim("EUR_candidate_lgd_allTargets.relateOut.txt", header=F, sep=";")
prob <- read.delim("EUR_proband_lgd_allTargets.relateOut.txt", header=F, sep=";")
sib <- read.delim("EUR_sibling_lgd_allTargets.relateOut.txt", header=F, sep=";")

cand$age <- (cand$V9 + cand$V10)/2
prob$age <- (prob$V9 + prob$V10)/2
sib$age <- (sib$V9 + sib$V10)/2

cand$type = "Candidate_LGD"
prob$type = "Proband_LGD"
sib$type = "Sibling_LGD"

dat <- data.frame(Variant=c(cand$type, prob$type, sib$type), Age=c(cand$age, prob$age, sib$age))
dat_1k <- dat[ dat$Age <1000,]
dat_100 <- dat[ dat$Age <100,]
dat_25 <- dat[ dat$Age <25,]
dat_10 <- dat[ dat$Age <10,]



my_comparisons <- list( c("Candidate_LGD", "Proband_LGD"), 
                        c("Proband_LGD", "Sibling_LGD"), c("Candidate_LGD", "Sibling_LGD") )

col_palette = c("#FC4E07", "#E7B800", "#00AFBB")

p <- ggboxplot(dat, x = "Variant", y = "Age", color = "Variant",
               palette = col_palette, add = "jitter")  + xlab("Variant (full)") +
    stat_compare_means(comparisons = my_comparisons) + stat_compare_means(label.y = 50000) 

p_1k <- ggboxplot(dat_1k, x = "Variant", y = "Age", color = "Variant",
               palette =col_palette, add = "jitter")  + xlab("Variant (age < 1,000 generation)") +
    stat_compare_means(comparisons = my_comparisons) + stat_compare_means(label.y = 500) 

p_100 <- ggboxplot(dat_100, x = "Variant", y = "Age", color = "Variant",
               palette =col_palette, add = "jitter")  + xlab("Variant (age < 100 generation)") +
    stat_compare_means(comparisons = my_comparisons) + stat_compare_means(label.y = 50) 

p_25 <- ggboxplot(dat_25, x = "Variant", y = "Age", color = "Variant",
               palette =col_palette, add = "jitter")  + xlab("Variant (age < 25 generation)") +
    stat_compare_means(comparisons = my_comparisons) + stat_compare_means(label.y = 15) 
p_10 <- ggboxplot(dat_10, x = "Variant", y = "Age", color = "Variant",
               palette =col_palette, add = "jitter")  + xlab("Variant (age < 10 generation)") +
    stat_compare_means(comparisons = my_comparisons) + stat_compare_means(label.y = 5) 


pdf("cmp_age_lgd.pdf", width=18, height=6)
grid.arrange(p, p_25, p_10, nrow=1)
dev.off()




