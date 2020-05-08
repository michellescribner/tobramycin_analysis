# chi squared goodness of fit test to determine if fusA1 snps reported in Scribner et all
# are distributed randomly across the length of fusA

# Are mutations distributed among domains differently than expected based on number of mutations and domain size?
domain <- read.csv(file="/Users/mrs/Documents/tob_paper_data/fusA1_mutation_distribution/fusA_dist_domain.csv", header=TRUE)
head(domain)
x <- domain$count
p <- domain$p
chisq.test(x, p = domain$p)

#Chi-squared test for given probabilities
#data:  x
#X-squared = 36.956, df = 4, p-value = 1.839e-07

# Are mutations distributed among 100 amino acid intervals differently than expected based on number of mutations?
b100 <- read.csv(file="/Users/mrs/Documents/tob_paper_data/fusA1_mutation_distribution/fusA_dist_100.csv", header=TRUE)
head(b100)
x <- b100$count
chisq.test(x, p = b100$p)

#Chi-squared test for given probabilities
#data:  x
#X-squared = 54.175, df = 7, p-value = 2.172e-09

# Are mutations distributed among 100 random bins differently than expected based on number of mutations and size of the bins?
rand <- read.csv(file="/Users/mrs/Documents/tob_paper_data/fusA1_mutation_distribution/fusA_dist_rand.csv", header=TRUE)
head(rand)
x <- rand$count
chisq.test(x, p = rand$p)

#Chi-squared test
#data:  x
#X-squared = 41.957, df = 4, p-value = 1.703e-08
