#6970 Assignment 1- Full R Script
#Sichong Xu

#Load the required packages
library(ggplot2)
library(car)
library(forecast)

#Problem 1#####

#Load in the design matrix X
design.matrixX <- read.csv("model.matrix.csv")
X <- as.matrix(design.matrixX)

#Question i a)
#M = X'X
M <- t(X)%*%X

#To find all eigenvalues and eigenvectors for matrix M
ev <- eigen(M)
eigenvalues <- ev$values
eigenvectors <- ev$vectors
eigenvalues
eigenvectors

#Question i b)
#To find an important property of matrix M.
solve(M) #M inverse
det(M)
#As the determinant of M is considered as 0, M is non-invertible, redundancy in matrix


#Question i c)
#If x is an eigenvector of M, then Mx = λx
egvec1 <- as.numeric(eigenvectors[,1])  #The first eigenvector of all
egv1 <- as.numeric(eigenvalues[1])#Its corresponding eigenvalue

left <- M%*%egvec1 #Left side of equation
right <- egv1%*%egvec1  #right side of equation

#Check if left side of the equation is equal to the right side
left
right
#Left and right side are the identical, so the the first eigenvector is indeed an eigenvector of M.


#Question i d)
#Use the rate data from "Puromycin" dataset, and attach to design matrix X.
data("Puromycin")
Xdata<- as.data.frame(cbind(Puromycin[2],X))

#Fit a linear regression model by regressing rate on conc, stateuntreated, conc.stateuntreated, statetreated.
frate4 <- lm(rate ~ conc + stateuntreated + conc:stateuntreated + statetreated, data = Xdata)
summary(frate4)
model.matrix(frate4)

#Check the estimated regression coefficients of the fitted model “frate4”.
coef(frate4)
#There is a missing value for "statetreated".

Xdata$stateuntreated
Xdata$statetreated
#Both "stateuntreated" and "statetreated" are binary variables. Variable "stateuntreated" has "1" represents the cell is not treated with Puromycin, "0" represents the treated status. On the contrary, "statetreated" represents the opposite situation.
#Therefore, only one of these binary variables is sufficient to indicate the status of treatment. To avoid the redundancy of fitted model, the regression coefficient of "statetreated" is treated as missing value "NA".

#A desired modification on matrix X is to remove the column of variable "stateuntreated", then a more plausible linear regression model can be built.
Xm<- Xdata[,-6]
fm<- lm(rate ~ conc + stateuntreated + conc:stateuntreated, data = Xm)
summary(fm)
coef(fm)


#Question ii) 
#Eye examination for any potential relationship.
plot(Puromycin)

#Generate a potential correlation plot.
ggplot(Puromycin, aes(x = conc, y = rate, shape = state, color = state)) + geom_point(cex=2) +
  ggtitle ("Rate vs. Concentration (Puromycin treated/untreated)")

#Fit a linear model, regressing rate on other two variables.
f0 <- lm(rate ~., data = Puromycin)
summary(f0)
#Base on the p-value, conc is a significant variable in this model.

par(mfrow=c(2,2))
plot(f0)
par(mfrow=c(1,1))
#Residual plots show non-linearity. 

ggplot(Puromycin, aes(x = conc, y = rate, shape = state, color = state)) + 
  geom_point(cex = 2) + 
  geom_smooth(method = lm, mapping = aes(y = predict(f0, Puromycin)))
#In the plot, we can observe this model has two parallel lines with same slope, which assumes that the effect of contraction is the same for treated and untreated cells. 
#As the assumption does not meet the previous observation of Puromycin dataset, it is clearly that data transformations and interaction effects should be introduced to produce a better model.

ggplot(Puromycin, aes(x = conc, y = log(rate), shape = state, color = state)) + geom_point(cex=2)
ggplot(Puromycin, aes(x = conc, y = sqrt(rate), shape = state, color = state)) + geom_point(cex=2)
#No linearity is found between "log(rate)" or "sqrt(rate)" and "conc".

#But still try with log-transformation on the response variable "rate" to build a model.
f1 <- lm(log(rate) ~., data = Puromycin)
summary(f1)

#Square root transformation on "rate".
f2 <- lm(sqrt(rate) ~., data = Puromycin)
summary(f2)

summary(f0)$adj; summary(f1)$adj; summary(f2)$adj
#Clearly, f0 with no modification on "rate" has the best adjusted R-squared value so far. 

ggplot(Puromycin, aes(x=log(conc), y=rate, shape=state, color=state)) + geom_point(cex=2)
ggplot(Puromycin, aes(x=sqrt(conc), y=rate, shape=state, color=state)) + geom_point(cex=2)
#Log transformation on the predictor "conc" successfully reveals linear relationship between "rate" and "log(conc)" with different status of "state".

#Fit model with log-transformation on the predictor (conc)
f3 <- lm(rate ~ log(conc) + state,  data = Puromycin)
summary(f3)
summary(f0)$adj; summary(f3)$adj
#Log-transformation on "conc" gives higher adjusted R-squared value, f3 is a better model.

#Try with square root transformation on the predictor (conc) for a new model.
f4 <- lm(rate ~ sqrt(conc) + state,  data = Puromycin)
summary(f4)
summary(f3)$adj; summary(f4)$adj
AIC(f3); AIC(f4)
#Model f3 has higher adjusted R-squared value.

#Add an interact effect "conc:state" to model f3.
f5 <- lm(rate ~ log(conc) + state + conc:state,  data = Puromycin)
summary(f5)
summary(f3)$adj; summary(f5)$adj
#Even f5 has a slightly higher adjusted R-squared value, the p-value associated with both interaction terms are not statistically significant.

#Try with interact effect "log(conc):state".
f6 <- lm(rate ~ log(conc) + state + log(conc):state,  data = Puromycin)
summary(f6)

summary(f3)$adj; summary(f5)$adj; summary(f6)$adj
AIC(f3); AIC(f5); AIC(f6)
#The model f6 would be considered as the most superior model after all kinds of transformations, which has the largest adjusted R-squared value, and the smallest Akaike information Criteria (AIC). 
#Also, the interaction term "log(conc):stateuntreated" is considered as significant. Therefore, the effect of log transferred concentration and the effect of treatment are dependent on one another.  

#The final best model can be designed.
final1 <- f6

#To extract the estimated model coefficients for identifying the equations.
round(coef(final1),2)
#Generate the equation for the final best regression model. 
#formula = rate ~ log(conc) + state + log(conc):state, data = Puromycin
#rate = 209.19 + 37.11 * log(conc) - 44.61 * [stateuntreated] - 10.13 * log(conc) * [stateuntreated]

#visualize the regression model.
ggplot(Puromycin, aes(x = log(conc), y = rate, shape = state, color = state)) + 
  geom_point(cex = 2) + 
  geom_smooth(method = lm, mapping = aes(y = predict(final1, Puromycin))) + 
  ggtitle ("Rate vs. log(Concentration) and Puromycin Treatment")
#Noticed that two lines have different slopes, which indicates that regarding the treatment state of cells, regression model f6 assumes different increase rates of the reaction rate as response to the predictor concentration (conc).

#Equations for  two straight lines in model.
#For state = untreated: rate = 209.19 + 37.11 * log(conc) - 44.61 * 1 - 10.13 * log(conc) * 1
#                       rate = 164.58 + 26.98 * log(conc)
#For state = treated: rate = 209.19 + 37.11 * log(conc) - 44.61 * 0 - 10.13 * log(conc) * 0
#                     rate = 209.19 + 37.11 * log(conc)


#Question iii) 
#Let's visualize the model to check its goodness-of-fit.
#First take a look at the residual plots of each model.
plot(f0, which = 1)
plot(f5, which = 1)
plot(final1, which = 1)
#The residuals of final1 is more randomly and more equally distributed around a horizontal dotted line than other models.
par(mfrow=c(2,2))
plot(final1, main = "Quality Check for the Final Best Model")
par(mfrow=c(1,1))


#Question iv) 
#To apply the best linear regression model fitted above to predict the reaction rate at at a concentration level of 0.15 ppm.
newdata <- data.frame(conc = c(0.15, 0.15), state = c("untreated", "treated"))
predict(final1, newdata)
predict(final1, newdata, interval = 'confidence')



#Problem 2#####
#Question i)
#Read in the input data.
data.train <- read.csv('dat.train.csv')

#Fit a simple linear regression model to test the hypothesis
res <- lm(weight ~ habit, data = data.train)
summary(res)
summary(res)$coeff[,4]

#Question ii)
#Check the missing value each variable.
apply(is.na(data.train), 2, sum)

#Since it is not ideal to exclude missing data or treat missing data as 0, replace NAs with the mean of "fage"
data.train$fage[is.na(data.train$fage)] <- mean(data.train$fage, na.rm = TRUE)

#Convert data type to factor
data.train$mature <- as.factor(data.train$mature) 
data.train$premie <- as.factor(data.train$premie)
data.train$marital <- as.factor(data.train$marital)
data.train$gender <- as.factor(data.train$gender)
data.train$habit <- as.factor(data.train$habit)
data.train$whitemom <- as.factor(data.train$whitemom)

#Plot the correlation between individual variables.
plot(data.train)
#Plot the bivariate relationships between the combination of numeric variables
scatterplotMatrix(~ fage + mage + weeks + visits + gained + weight,
                  data = data.train, main = "Scatterplot Matrix",
                  smooth = list(spread = T, col.smooth = "black", col.spread = "black", 
                              lty.smooth = 2, lwd.smooth = 3, lty.spread = 3, lwd.spread = 2),
                  col = c("red"),)
#Both visualizations show certain linear relationship between "wright" and "weeks".

#Applying log transformation to "weight", to see if more stable correlation can be determined.
scatterplotMatrix(~ fage + mage + weeks + visits + gained + log(weight),
                  data = data.train, main = "Scatterplot Matrix with log(weight)",
                  smooth = list(spread = T, col.smooth = "black", col.spread = "black", 
                              lty.smooth = 2, lwd.smooth = 3, lty.spread = 3,lwd.spread = 2),
                  col = c("red"),)

#Fit the first model
ffull <- lm(weight ~., data = data.train) 
summary(ffull)
#By looking at the p-value, many variables is not significant.

pvalues <- summary(ffull)$coeff[,4]
vnames <- c(rownames(summary(ffull)$coeff)[pvalues <.1])
vnames
#The variable has small p-value ([pvals <.1])) is determined for following model building.

#Automatic model selection.
select_back <- step(ffull, direction = c("backward"))  
select_forw <- step(ffull,direction = c("forward"))
select_stepwise <- step(ffull, direction = c("both")) 
summary(select_back)
summary(select_forw)
summary(select_stepwise)
AIC(select_back);AIC(select_forw); AIC(select_stepwise)
#Stepwise selection has the best result.

#Stepwise regression is performed using AIC to select variables for model.
fsimple <- step(ffull, direction = "both") 
summary(fsimple)

#Fit a model by regressing log(weight) to other variables, and determine the simplified model with stepwise selection.
ffull_log <- lm(log(weight) ~., data = data.train) 
summary(ffull_log)
fsimple_log <- step(ffull_log, direction = "both") 
summary(fsimple_log)

#Compare models.
summary(ffull)$adj; summary(fsimple)$adj; summary(ffull_log)$adj; summary(fsimple_log)$adj
AIC(ffull);AIC(fsimple); AIC(ffull_log); AIC(fsimple_log)
#fsimple_log is the best one among 4 models

#Plot residuals to exam the quality of model.
plot(fsimple_log, which = 1)
crPlots(fsimple_log)
#Non-linearity is found between "log(weight)" and "weeks" in the residuals plot.

#Log transformation on "weeks" is applied to fit a better model.
ffull_logweeks <- lm(log(weight) ~ fage + mage + mature + 
                       log(weeks) + premie + visits + marital + 
                       gained + gender + habit + whitemom , data = data.train)
summary(ffull_logweeks)
#Stepwise selection (AIC).
fsimple_logw <- step(ffull_logweeks, direction = "both")
summary(fsimple_logw)
#Plot residuals.
plot(fsimple_logw, which = 1)
crPlots(fsimple_logw)
#Plots shows the log transformation on "weeks" is not sufficient enough to fix the non-linearity.

#A log-transferred "weeks" is add to the above model.
ffull_logweeks2 <- lm(log(weight) ~. + log(weeks), data = data.train) 
summary(ffull_logweeks2)

fsimple_logw2 <- step(ffull_logweeks2, direction = "both") 
summary(fsimple_logw2)

plot(fsimple_logw2, which = 1)
crPlots(fsimple_logw2)
#Great improvement can be seen in component+residual plots.

summary(fsimple_log)$adj; summary(fsimple_logw)$adj; summary(fsimple_logw2)$adj
AIC(fsimple_log); AIC(fsimple_logw); AIC(fsimple_logw2)
#fsimple_logw2 is the best model so far.

#Try with adding polynomial effect to "weeks" with "log(weight)".
fplm <- lm(log(weight) ~ fage + mage + mature + 
             poly(weeks, 2) + premie + visits + marital + 
             gained + gender + habit + whitemom, data = data.train)
fsimple_plm<- step(fplm, direction = "both") 
summary(fplm)
summary(fsimple_plm)

fplm4 <- lm(log(weight) ~ fage + mage + mature + 
             poly(weeks, 4) + premie + visits + marital + 
             gained + gender + habit + whitemom, data = data.train)
fsimple_plm4<- step(fplm4, direction = "both")
summary(fsimple_plm4)

fplm5 <- lm(log(weight) ~ fage + mage + mature + 
              poly(weeks, 5) + premie + visits + marital + 
              gained + gender + habit + whitemom, data = data.train)
fsimple_plm5<- step(fplm5, direction = "both")
summary(fsimple_plm5)

AIC(fsimple_logw2); AIC(fsimple_plm4); AIC(fsimple_plm5)
#Fourth power on weeks gives the best model

#Remove the non-significant variable "premie".
ftest <- lm(log(weight) ~ poly(weeks, 4) + marital + 
              gained + gender + habit, data = data.train)
summary(ftest)
AIC(fsimple_plm4); AIC(ftest)
#Result shows keeping "premie" as the covariable, model has slightly better adjusted R-squared value and AIC.

#Use Anova to compare the models.
anova(fsimple_log, fsimple_plm4)

anova(fsimple_logw2, fsimple_plm4)
#Significant difference between two models can be observed based on the p-value in ANOVA test, it is meaningful to adopt the polynomial effect poly(weeks, 4).
#The current model with log transformation on "weight" and fourth power on "weeks" is the best model.

#The final/ best multiple linear regression model can be determined.
final <- fsimple_plm4

#Quality check with residuals plots.
par(mfrow=c(2,2))
plot(final)
par(mfrow=c(1,1))

crPlots(final)

#Whether "habit" is an important predictor in the final model.
summary(final)$coeff[,4]
#"habit" is an important predictor in this model, the conclusion id i) does not agree with the final model,


#Question iii) a.
#Load the data.
data.test <- read.csv('dat.test.csv')

#Predict birthweight for observations in the test dataset with the final model.
pred <- exp(predict(final, data.test))
pred

pred_interval <- exp(predict(final, data.test, interval="prediction"))

data.test1 <- data.test[,-9]
data.test1$predweight <- exp(predict(final, data.test1))
data.test1$predweight

#Question iii) b.
#Predict birthweight with other model fitted above.
pre_logw2 <- exp(predict(fsimple_logw2, data.test))

#Compare the predicted value with the true value obtain from test dataset
true <- data.test$weight

#Measure the accuracy of 2 prediction models, calculate the RMSE: Root Mean Squared Error.
accuracy(true,pre_logw2)
accuracy(true, pred) 

#Compute the Mean Squared Prediction Error
mse_pre_logw2 <- mean((true - pre_logw2)^2)
mse_pred <- mean((true - pred)^2)
mse_pre_logw2 
mse_pred 

#Visualizing the comparison.
plot(true, pred) #True vs. predicted value
abline(0,1,col='red') #Add a line y=x

#Draw a better scatter plot with ggplot2.
res <- data.frame(prediction = pred, true = true)
ggplot(res, aes(true, prediction)) + geom_point() + geom_abline(intercept = 0, slope = 1,col = 'red') +
  ggtitle ("Comparison Between Observed and Predicted Birthweight")



