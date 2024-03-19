library(forecast)
library(quantmod)
library(TSA)
library(tseries)  
library(urca)   ## aquí esta el test de dickey-fuller
#library(lmtest)

imae <- imae_variacion_acumulada$IMAE
plot(imae)
imae_inter <- imae_variacion_interanual$IMAE
imae_ciclo <- imae_medidas_absolutas$C38


#convirtiendo los datos en un objeto time series
imae_ts <- ts(imae, start = c(2001, 01), frequency = 12)

#Grafica del IMAE
plot(imae_ts, main = "Serie Variación Acumulada - IMAE", lwd = 1,
     xlab = "Años", ylab = "IMAE")


# prueba de raiz unitaria Dickey-Fuller

# H_0: No hay estacionariedad || Hay raíz unitaria
#H_a: Hay estacionariedad     || No hay raiz unitaria
# Rechazamos H_0 si p-value < alpha o t < v. critic

#Si la serie presenta tendencia, incluir tendencia e intercepto
#si no tiene tendencia y su media no es cero, incluir solo intercepto
#si fluctua en torno a su media (media cero), no incluir nada  

# las siguientes pruebas son ambas con constante e intercepto
adf.test(imae_ts, k=0) #D-F estandar  (según el test, no hay estacionariedad)
df1 <- ur.df(imae_ts, type = "trend", lags = 0)
summary(df1)
# prueba solo con intercepto
df2 <- ur.df(imae_ts, type = "drift", lags = 0) 
summary(df2) # segun el test, no hay estacionariedad

#prueba sin incluir intercepto y tendencia
df3 <- ur.df(imae_ts, type = "none", lags = 0)
summary(df3) #segun el test, no hay estacionariedad


# Descomosición de la serie para observar sus componentes (tendencia, estacionalidad...)
plot(decompose(imae_ts))

#Pruebas de autocorrelación para los datos de la serie
acf(imae_ts, type = "correlation", main = "Autocorrelación de los datos")
acf(imae_ts, type = "partial", main = "Autocorrelación parcial de los datos")

# otras funciones para pruebas de autocorrelación
acf(imae_ts, main = "ACF de las observaciones", lag.max = 50)
pacf(imae_ts, main = "PACF de las observaciones", lag.max = 50)

#numero de veces que necesitamos diferenciar la serie para volverla estacionaria
ndiffs(imae_ts)
imae_diff <- diff(imae_ts)
plot(imae_diff)

#Prueba de Dickey fuller a la serie diferenciada
adf.test(imae_diff,k = 0)
diff_df <- ur.df(imae_diff, type = "none", lags = 0)
summary(diff_df) # la serie diferenciada es estacionaria

#acf a la serie differenciada
acf(imae_diff, main = "ACF de serie diferenciada", lag.max = 50)
pacf(imae_diff, main = "PACF de serie diferenciada", lag.max = 50)


# crear modelo arima 1 0 1 :
arima101 <- arima(imae_ts, order = c(1,0,1))
arima101

#test de coeficientes

coeftest(arima101)

#arima 1 1 1:
arima111 <- arima(imae_ts, order = c(1,1,1))
arima111

coeftest(arima111)

#arima 2 1 2:
arima212 <- arima(imae_ts, order = c(2,1,2))
arima212

coeftest(arima212)

##### predicción
pred212 <- predict(arima212, n.ahead = 12)
pred212
pred101 <- predict(arima101, n.ahead = 12)
pred101
pred111 <- predict(arima111, n.ahead = 12)
pred111

#comparando serie original con la de los modelos
plot(imae_ts, 
     main = "Serie Original vs Serie ARIMA",
     xlab = "Años", ylab = "IMAE", lwd = 1)
lines(fitted.values(arima212), col = "red", lwd = 1)
lines(pred212$pred, col= "red", lwd = 3)
lines(fitted.values(arima111), col = "blue", lwd = 1)
lines(pred111$pred, col= "blue", lwd = 3)
lines(fitted.values(arima101), col = "green", lwd = 1)
lines(pred101$pred, col = "green", lwd = 3)

plot(arima212)


resid212 <- resid(arima212)
plot(resid212)


###################################################################
######### probando el modelo arima(2,1,2)

# Se tomarán los datos hasta diciembre 2022
imae_truncado <- ts(imae, start = c(2001,1), end = c(2022,12), frequency = 12)
plot(imae_truncado)

arima212_truncado <- arima(imae_truncado, order = c(2,1,2))
arima212_truncado
coeftest(arima212_truncado)


### Residuos del modelo arima
resid212 <- resid(arima212_truncado)
plot(resid212)

### observar la autocorrelación de los residuos
#   y comprobar que haya independencia entre ellos,
#   esto muestra mayor exactitud del modelo
acf(resid212, main = "ACF de los residuos - ARIMA(2,1,2)", 
    lag.max = 150)

### Prueba Box-Ljung para descartar autocorrelación
#   en los residuos
### H_0: p_0 = p_1 = ... = p_n = 0
#   H_a: p_i != 0 para algún i 
#   rechazamos H_0 si p < 0.05
Box.test(resid212, lag = 6, type = "Ljung") # m ~ Ln(264), Q(6)
#segun la prueba realizada, no hay dependencia entre los resid

#pronóstico para 24 meses == 2 años
pred212_truncado <- predict(arima212_truncado, n.ahead = 24)
pred212_truncado

#construyendo los intervalos de confianza:
v.critic <- qnorm(0.975) # valor critico al 95%
lim_inf <- pred212_truncado$pred - v.critic*pred212_truncado$se
lim_sup <- pred212_truncado$pred + v.critic*pred212_truncado$se

pronostico <- data.frame("fecha" = time(pred212_truncado$pred),
                         "predicción" = pred212_truncado$pred,
                         "lim inf" = lim_inf,
                         "lim sup" = lim_sup)
pronostico

#graficando serie original vs modelo
plot(imae_ts, lwd = 1)
lines(imae_truncado, col = "blue", lwd = 2)
lines(fitted.values(arima212_truncado), col = "red", lwd = 1)
lines(pred212_truncado$pred, col= "red", lwd = 1)
lines(lim_inf, col = "brown", lwd = 2)
lines(lim_sup, col = "brown", lwd = 2)

#Vamos a calcular el error de mis datos estimados de 2023
# con los reales, es decir, error = | x^ - x |:
data23 <- imae[265:275]
data23_estimada <- pred212_truncado$pred[1:11]
error212 <- abs(data23_estimada - data23)
error212
mean(error212)
#graficando los errores:
plot(error212)

#Calculando el MSE (error cuadrático medio):
data <- data.frame("observado" = imae[1:264], "estimado" = fitted.values(arima212_truncado))
mse212 <- mean((data$estimado - data$observado)^2)
mse212

# el MSE también se puede obtener directo de los residuos del modelo
mse212 <- mean(arima212_truncado$residuals^2)
mse212

###################################################################
######### probando el modelo arima(1,1,1)

#Se tomarán los datos hasta diciembre 2022
imae_truncado <- ts(imae, start = c(2001,1), end = c(2022,12), frequency = 12)
plot(imae_truncado)

arima111_truncado <- arima(imae_truncado, order = c(1,1,1))
arima111_truncado
coeftest(arima111_truncado)
pred111_truncado <- predict(arima111_truncado, n.ahead = 12)
pred111_truncado

plot(imae_ts, lwd = 1)
lines(imae_truncado, col = "blue", lwd = 2)
lines(fitted.values(arima111_truncado), col = "red", lwd = 1)
lines(pred111_truncado$pred, col= "red", lwd = 1)


####################################################################
################ usando funcion auto.arima

#buscando un modelo sarima
autosarima <- auto.arima(imae_truncado)
autosarima
coeftest(autosarima)

pred_sarima <- predict(autosarima, n.ahead = 24)
pred_sarima

plot(imae_ts, lwd = 1)
lines(imae_truncado, col = "blue", lwd = 2)
lines(fitted.values(autosarima), col = "red", lwd = 1)
lines(pred_sarima$pred, col= "red", lwd = 1)

#buscando un modelo arima
autoarima <- auto.arima(imae_truncado, seasonal = FALSE)
autoarima
coeftest(autoarima)

pred110 <- predict(autoarima, n.ahead = 24)
pred110

plot(imae_ts, lwd = 1)
lines(imae_truncado, col = "blue", lwd = 2)
lines(fitted.values(autoarima), col = "red", lwd = 1)
lines(pred110$pred, col= "red", lwd = 1)
