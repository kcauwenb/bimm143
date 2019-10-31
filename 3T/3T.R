# 3T: Data Visualization

x<-rnorm(1000)
mean(x)
sd(x)
summary(x)
boxplot(x)
hist(x)
# adds tassels at bottom of histogram, shows where data specifically is
rug(x)

# break

# Section 2: Customizing plots
# 2A. Line plot
weight_chart<-read.table("bimm143_05_rstats/weight_chart.txt", header = TRUE)
plot(weight_chart$Age, weight_chart$Weight,type="o",pch=15,cex=1.5,lwd=2,
     ylim=c(2,10),xlab="Age (months)",ylab="Weight (kg)",
     main="Overplotted points/lines: Age vs. Weight",col="blue")

# 2B. Barplot
mouse<-read.table("bimm143_05_rstats/feature_counts.txt",header=TRUE,sep = "\t")
# equivalent to above
mouse<-read.delim("bimm143_05_rstats/feature_counts.txt")
# set margins
par(mar=c(3.1,11.1,4.1,2))
# las=1 sets axis labels horizontal, look at ?par
barplot(mouse$Count, horiz = T, names.arg=mouse$Feature,
        main="Mouse Genomic Features",las=1,xlim=c(0,80000))

# 2C. Histograms
x <- c(rnorm(10000),rnorm(10000)+4)
hist(x,breaks = 80)


# Section 3: Using color in plots
# 3A. Providing color vectors
male_female<-read.delim("bimm143_05_rstats/male_female_counts.txt")
par(mar=c(10,5,5,5))
barplot(male_female$Count, names.arg = male_female$Sample,
        las=2, col=rainbow(nrow(male_female)))
# col argument allows alternationg colors
barplot(male_female$Count, names.arg = male_female$Sample,
        las=2, col=c("blue","red"))

# 3B. Coloring by value
genes<-read.delim("bimm143_05_rstats/up_down_expression.txt")
# match state categories with colors
levels(genes$State)
palette(c("blue","grey","red"))
# plot c1 against c2, with colors being determined by state
plot(genes$Condition1,genes$Condition2,col=genes$State,
     xlab="Condition 1", ylab="Condition 2")

# 3C. Dynamic use of color
meth<-read.delim("bimm143_05_rstats/expression_methylation.txt")
# create vector of blue hues according to density
col_density=densCols(meth$gene.meth,meth$expression)
# plot data, pch=20 changes unfilled circle to filled
plot(meth$gene.meth,meth$expression,col=col_density,
     pch=20)
# now, restrict to genes with higher expression levels
hi_exp_indices<-meth$expression>0
plot(meth$gene.meth[hi_exp_indices], 
     meth$expression[hi_exp_indices])
dcols <- densCols(meth$gene.meth[hi_exp_indices],
                  meth$expression[hi_exp_indices])
plot(meth$gene.meth[hi_exp_indices], 
     meth$expression[hi_exp_indices], 
     col = dcols, pch = 20)
# now, have the colors go between blue, green, red and yellow
dcols.custom <- densCols(meth$gene.meth[hi_exp_indices], 
                         meth$expression[hi_exp_indices],
                         colramp = 
                        colorRampPalette(c("blue","green","red","yellow")))

plot(meth$gene.meth[hi_exp_indices], 
     meth$expression[hi_exp_indices], 
     col = dcols.custom, pch = 20)
