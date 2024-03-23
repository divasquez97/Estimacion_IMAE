library(forecast)
library(ggplot2)
library(TSA)
library(quantmod)
library(tseries)
library(stats)
library(urca)
library(lmtest)
library(dynlm)
library(timeSeries)
library(FinTS)
library(rugarch)
library(fGarch)
library(rmgarch)

imae <- imae_variacion_acumulada$IMAE
plot(imae)

#tomando datos hasta diciembre 2022
#convirtiendo los datos en un objeto time series
imae_ts <- ts(imae, start = c(2001, 01), frequency = 12)
imae_truncado <- ts(imae, star = c(2001,01), end = c(2022,12), frequency = 12)

#arima 2 1 2:
arima212 <- arima(imae_truncado, order = c(2,1,2))
arima212

coeftest(arima212)

pred212 <- predict(arima212, n.ahead = 24)
pred212

#comparando serie original con la del modelo
plot(imae_ts, 
     main = "Serie Original vs Serie ARIMA",
     xlab = "Años", ylab = "IMAE")
lines(imae_truncado, col = "blue", lwd = 2)
lines(fitted.values(arima212), col = "red", lwd = 2)
lines(pred212$pred, col = "brown", type = "l", lwd = 2)

################################################################
####################### ARCH - GARCH ###########################

#residuos al cuadrado del modelo arima(212)
rescuad212 <- resid(arima212)^2
# chartSeries(rescuad212, 
#             theme = chartTheme("white"))
plot(rescuad212, main = "Residuos^2 ARIMA(2,1,2)", 
     xlab = "Años", ylab = "Residuos", col = "blue4")



#revisando la autocorrelacion de los residuos
acf(rescuad212, main = "ACF residuos^2", lag.max = 100)
pacf(rescuad212, main = "PACF residuos^2", lag.max = 100)

#Prueba de efecto ARCH con un rezago para mis datos IMAE
ArchTest(imae_ts, lags = 1, demean = TRUE) # la prueba dice que hay efecto arch

#especificar el modelo garch deseado
# en este caso se quiere un ARIMA(2,1,2) + GARCH(1,1) 
mlgarch <- ugarchspec(mean.model = list(armaOrder = c(2,2), 
                                        arfima = T), 
                      fixed.pars = list(arfima = 1))
mlgarch

#estimar el modelo garch
fit_garch <- ugarchfit(spec = mlgarch, data = imae_truncado)
fit_garch@fit$coef

#pronosticos del modelo para 24 meses == 2 años
garch_forecast <- ugarchforecast(fit_garch, n.ahead = 11)
garch_forecast

plot(garch_forecast)

#Calcular los errores, error = |x^ - x|
data23 <- imae[265:275]
data23_estimada <- garch_forecast@forecast$seriesFor[1:11]
error_garch <- abs(data23_estimada - data23)
error_garch

#graficando los errores:
plot(error_garch)

#MSE
residuos_garch <- garch_forecast@model$modeldata$residuals
mse_garch <- mean(residuos_garch^2)
mse_garch

#Guardando los valores estimados del modelo para poder
# graficarlos de forma adecuada
fittedvalues_garch <- ts(fit_garch@fit$fitted.values, 
                         start = c(2001,01), frequency = 12)
#lo mismo para los valores pronosticados
pronost_garch <- ts(garch_forecast@forecast$seriesFor, 
                    start = c(2023,01), frequency = 12)

#lo anterior solamente lo necesito para graficar
# la serie original vs la del modelo bien indexada 

#Graficar el modelo vs serie original
plot(imae_ts, main = "Serie original vs GARCH", 
     type = "l", lwd = 2)
lines(fittedvalues_garch, col = "red", lwd = 2)
lines(pronost_garch, col = "red", lwd = 2)


################################################################
################################################################
####################### ARCH - GARCH ###########################
chartSeries(imae_ts, name = "IMAE - V.A")

#Diferenciar la serie original para luego aplicar un ARMA(2,2)
imae_diff <- diff(imae_ts)
plot(imae_diff, main = "Serie diferenciada" , lwd = 2)

#modelo arma(2,2)
arma22 <- arima(imae_diff, order = c(2,0,2))
arma22
coeftest(arma22)

#graficando el modelo 
plot(imae_diff, lwd=2)
lines(fitted.values(arma22), lwd = 2, col = "red")

#calcular residuos al cuadrado con el arma(22) de mi serie diff
rescuad <- resid(arma22)^2
rescuad
chartSeries(rescuad)

##### lo anterior lo he hecho debido a que el modelo garch
# necesita asociar un modelo arma, pero el modelo que se ajusta
# mejor a mis datos es un arima(2,1,2), entonces mi idea es
# diferenciar 1ero mi serie y luego aplicar un arma(2,2)
# que será el modelo que asociaré a mi garch


#El efecto arch se refiere a que la varianza es heterocedastica
# y depende de los residuos^2 rezagados. Asocia un proceso AR

#En el efecto garch la varianza heterocedastica depende de
# los residuos^2 rezagados + varianza rezagada. Asocia un ARMA


#revisando la autocorrelacion de los residuos
acf(rescuad, main = "ACF residuos^2", lag.max = 100)
pacf(rescuad, main = "PACF residuos^2", lag.max = 100)

# prueba con los datos diferenciados
ArchTest(imae_diff, lags = 1, demean = TRUE) # no hay efecto arch


#especificar el modelo garch deseado
mlgarch <- ugarchspec(mean.model = list(armaOrder = c(2,2)))
mlgarch

#estimar el modelo garch
fit_garch <- ugarchfit(spec = mlgarch, data = imae_ts)
fit_garch

#estimar el garch a los datos diferenciados
fit_garch <- ugarchfit(spec = mlgarch, data = imae_diff)
fit_garch

#coeficientes del ARMA(2,2) + GARCH(1,1)
fit_garch@fit$coef

#varianza
garch_var <- fit_garch@fit$var
garch_var

#residuos
garch_residuos <- (fit_garch@fit$residuals)^2
plot(garch_residuos, type = "l")
lines(garch_var, col = "green")

#pronosticos del modelo
garch_forecast <- ugarchforecast(fit_garch, n.head = 12)
garch_forecast


plot(fit_garch)
plot(garch_forecast)

# Graficando el garch
plot(fit_garch@fit$fitted.values, main = "Serie original vs GARCH", 
     type = "l", lwd = 2, col = "red")
lines(imae, type = "l")
lines(garch_forecast@forecast)
lines(garch_forecast@forecast$seriesFor, col = "red")





