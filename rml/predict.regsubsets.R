# this function is used for predict regsubsets in leaps package.

predict.regsubsets = function(object, newdata, id, ...) {
        form  <-  as.formula(~.)
        mat  <-  model.matrix(form, newdata)
        coefi  <-  coef(object, id)
        xvars  <-  names(coefi)
        mat[, xvars] %*% coefi
}