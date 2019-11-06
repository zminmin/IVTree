# Instrumental Variable Tree (IVTree) Introduction

<!-- The _causalTree_ function builds a regression model and returns an _rpart_ object, which is the object derived from _rpart_ package, implementing many ideas in the CART (Classification and Regression Trees), written by Breiman, Friedman, Olshen and Stone. Like _rpart_, _causalTree_ builds a binary regression tree model in two stages, but focuses on estimating heterogeneous causal effect. -->

The _honest.IVTree_ function builds a honest instrumental variable tree model and returns an _rpart_ object. The function inherits many ideas in the CART (Classficiation and Regression Trees by Breiman, Friedman, Olshen and Stone) and the causalTree (causal tree by Athey and Imbens). By specifying the instrumental variable in the model, tree built by _IVTree_ corrects for potential endogeneity issues in observational data.


To install this package in R, run the following commands:

```R
install.packages("devtools")
library(devtools)
install_github("gweehwa/IVTree")
```

Example usage:

```R
library(IVTree)
tree <- honest.IVTree(y~ x1 + x2 + x3 + x4, data = data_tr, treatment = data_tr$T1, 
	treatment1 = $T, IV = data_tr$IV, split.Rule = "CT", 
	cv.option = "CT", split.Honest = T, cv.Honest = T, 
	split.Bucket = F, xval = 5, cp = 0, minsize = 20, propensity = 0.5)
                  
opcp <- tree$cptable[,1][which.min(tree$cptable[,4])]

opfit <- prune(tree, opcp)

rpart.plot(opfit)

```

For More details, please check out briefintro.pdf.

#### References
Guihua Wang, Jun Li, Wallace J. Hopp. <b>An Instrumental Variable Tree Approach for Detecting Heterogeneous Treatment Effects in Observational Studies.</b> [<a href="https://poseidon01.ssrn.com/delivery.php?ID=723089121105084082064005067083015023018031035064008038064102018004098091117000067094037026034111123061001010127099119014085106105082056047035121068113017124123069066040043085098007072027126103092105068089003088095011102077025105123083026069068065102004&EXT=pdf">link</a>]
