# fishers exact test to determine if mutations in fusA1 reported in Scribner et al occur more frequently in positions with identical amino acid identity across the five species studied than expected by chance.

######### fusA1 ##########
#             mutated	not mutated	
#identical    	31	      291
#not identical	15	      379

M <- as.table(rbind(c(31, 291), c(15, 379)))
dimnames(M) <- list(conservation = c("identical", "not identical"), 
                    mutation = c("mutated", "not mutated"))


(fusA_fisher <- fisher.test(M))
#Fisher's Exact Test for Count Data
#data:  M
#p-value = 0.001964
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 1.377815 5.464190
#sample estimates:
#odds ratio 
#  2.687957 
