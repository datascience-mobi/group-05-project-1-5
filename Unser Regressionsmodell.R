#library(pscl)


Reduce all genes to them of which we have expression, copynumber and knockdown data.
```{r 4_genes}
regression.genes = Reduce(intersect, list(rownames(expression), rownames(copynumber), rownames(kd.ceres)))[1:1000]
```

Create a dataframe with expression, copynumber, knockdown and mutation data from a HNC sample.
```{r 4_dataframe}
regression.genes = Reduce(intersect, list(rownames(expression), rownames(copynumber), rownames(kd.ceres)))

regression.mutation = sapply(1:length(regression.genes), function(a){
  ifelse(regression.genes[a] %in% mutation.HNC[[1]]$Hugo_Symbol, TRUE, FALSE)
})

regression.data = data.frame("genes" = regression.genes,
                             "expression" = expression.HNC[regression.genes, 1],
                             "copynumber" = copynumber.HNC[regression.genes, 1],
                             "kockdown" = kd.ceres.HNC[regression.genes, 1],
                             "mutation" = regression.mutation)

rm(regression.genes, regression.mutation)
```

## Train and test regression model
```{r 4_model, eval = FALSE, include = FALSE}
set.seed(123)
split = sample.split(regression.data$mutation, SplitRatio = 0.8)
training.set.rg = subset(regression.data, split == TRUE)
testing.set.rg = subset(regression.data, split == FALSE)
training.set.rg = regression.data1
testing.set.rg = regression.data2

train.model = glm(mutation~., family = binomial(link="logit"), data=training.set.rg)
summary(train.model)

anova(train.model, test="Chisq")

pR2(train.model)

model.results = predict(train.model, newdata = subset(testing.set.rg, select=c("genes", "expression","copynumber","kockdown")), type = "response") 
model.results = ifelse(model.results > 0.5,1,0)
```

## Compare different variable input
```{r 4_compare, eval = FALSE, include = FALSE}
model.results.knockdown = predict(train.model, newdata = subset(testing.set.rg, select="kockdown"), type = "response") # Hier wollte ich noch fragen welche Kombo ausprobiert werden soll
model.results.knockdown = ifelse(model.results.knockdown > 0.5,1,0)
```