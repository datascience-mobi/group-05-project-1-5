# K-means überlegungen

mean.cna.reference = apply(copynumber.reference, 1, mean) # generate mean copynumber values of the reference samples
cna.HNC.norm = as.data.frame(apply(copynumber.HNC, 2, function(x) {
  x - mean.cna.reference # use these mean values to normalise copynumber values from HNC samples
}))

expr.cna.genes = Reduce(intersect, list(rownames(expr.HNC.norm), rownames(cna.HNC.norm))) # choose only genes which we have copynumber and expression data from
km.expr = expr.HNC.norm[expr.cna.genes,]
km.cna = cna.HNC.norm[expr.cna.genes,]

dim(km.expr) == dim(km.cna) # Check if dimensions are equal

# Weg 1: x = expression, y = copynumber
km.data = lapply(1:length(ID), function(a) {
  out = data.frame(expression = km.expr[,a], copynumber = km.cna[,a])
  rownames(out) = rownames(km.expr)
  return(out)
})
names(km.data)[1:length(ID)] = ID

summary(kmeans(x = km.data[[1]], centers = 4))

plot(km.data[[2]]$expression, km.data[[2]]$copynumber, cex = 0.1)

rm(km.expr, km.cna, expr.cna.genes, km.data)

# Weg 2: k-means von expression und copynumber getrennt (ich weiß nicht ob das so viel Sinn macht ^^)
wss.exp = sapply(2:7, function(k) {
  kmeans(x = km.expr, centers = k)$tot.withinss # calculate wss for diferrent k centers for expression
})
plot(2:7, wss.exp, type = "b", pch = 19, xlab = "Number of clusters K (expression)", ylab = "Total within-clusters sum of squares")
# -> k = 3

wss.cna = sapply(2:7, function(k) {
  kmeans(x = km.cna, centers = k)$tot.withinss # calculate wss for diferrent k centers for copynumber
})
plot(2:7, wss.cna, type = "b", pch = 19, xlab = "Number of clusters K (Copynumber)", ylab = "Total within-clusters sum of squares")
# -> k = 5

km.expr.data = kmeans(x = km.expr, centers = 3, nstart = 10)
km.cna.data = kmeans(x = km.cna, centers = 4, nstart = 10)


summary(t(km.expr.data[["centers"]])) # center 1
summary(t(km.cna.data[["centers"]])) # center 2

km.expr.1 = names(km.expr.data$cluster[which(km.expr.data$cluster == 1)])
km.cna.2 = names(km.cna.data$cluster[which(km.cna.data$cluster == 2)])

high.exp.cna = Reduce(intersect, list(km.expr.1, km.cna.2))

rm(km.expr, km.cna, expr.cna.genes, wss.cna, wss.exp, km.cna.data, km.expr.data, km.expr.1, km.cna.2, high.exp.cna)
