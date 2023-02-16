pp_plot <- function(X, ag, bg, confidence, dist){
  #Cuantiles en el eje X y las empiricas en el Y
  Y = pnorm(sort(X), ag, bg) 
  X = seq(1/(n+1), 1, length.out = n)
  # puntos de la muestra
  plot(X,Y,
       main = "Gráfica PP de los residuos",
       xlab = "Probabilidades teóricas",
       ylab = "Probabilidades empíricas",
       pch = 19,
       cex = 0.5)
  # identidad
  abline(a = 0, b = 1, col = "red", lwd = 2)
  # bandas de confianza
  points(X,qbeta((1 - confidence)/2, 1:n, n + 1 - 1:n), 
         type = "l",
         lty = 2)
  points(X,qbeta((1 + confidence)/2, 1:n, n + 1 - 1:n), 
         type = "l",
         lty = 2)
}



plotRelative<-function(l, aG, bG, xL, xR, yL, yR, levels, n, xlab="", 
                       ylab="", main="Contornos de verosimilitud relativa"){
  x_vec = seq(from = xL, to = xR, length.out = n)
  y_vec = seq(from = yL, to = yR, length.out = n)
  R<-function(a, b){
    return(exp(l(a, b) - l(aG, bG)))
  }
  Rmat = matrix(nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      Rmat[i, j] = R(x_vec[i], y_vec[j])
    }
  }
  contour(x_vec,y_vec,Rmat,level=levels,xlab=xlab,ylab=ylab,main=main)
}


# Esta función recibe los vectores X, Y, n es el tamaño de la muestra y
# D es el determinante de K'K, donde K es la matriz de diseño.
# La función calcula la estimación de beta 1
B1 = function(X,Y,n,D){
  Prod = t(X) %*% Y
  return((n*Prod[1] - sum(X) * sum(Y)  ) / D)
}


# Esta función recibe los vectores X, Y, n es el tamaño de la muestra y
# b1 es la estimación de beta 1.
# La función calcula la estimación de beta 2
B0 = function(X,Y,n,b1){
  return((sum(Y)/n ) - (b1 * sum(X)/n))
}


# Esta función recibe el vector Y, la matriz de diseño K,
# n es el tamaño de la muestra, b1 es la estimación de beta 1 y 
# b2 es la estimación de beta 2.
# La función calcula la varianza estimada
vari = function(K,Y,n,b1,b2){
  B = matrix(c(b1,b2),nrow = 2, ncol=1)
  aux = ( t(Y - K %*% B)  %*%   (Y - K %*% B) ) 
  return( aux[1]  /    n)
}


#----------------------------------------------------------------------







Rpm = function(mu, n, t1, t2){
  return((t2 - (t1^2 / n))^(n/2)  / (t2 - 2*mu* t1 + n*mu^2 )^( n / 2))
}

Rpm1 = function(mu){
  return(Rpm(mu,n,t1,t2))
}


Rps = function(sigma,n,t1,t2){
  return( ((n*t2 - t1**2)**(n/2) * exp( (- (t2 - (t1**2 /n))/( 2* sigma**2) ) 
                                        + (n/2) ) )/ ((n*sigma)**n))
}


logver <- function(mu, sigma){
  return(-n*log(sigma) - (1/(2*sigma*sigma))*(n*mu*mu + t2 - 2*mu*t1))
}

logverx <- function(mu, sigma, t1, t2){
  return(-n*log(sigma) - (1/(2*sigma*sigma))*(n*mu*mu + t2 - 2*mu*t1))
}

AIC <- function(mg,sg,n,t1,t2){
  return(n* log( (2*3.1416 )) -2 * logverx(mg,sg,t1,t2)  + 4 )
}

AIC2 <- function(mg,sg,n,q1,q2){
  return(2*q1 + 2*n*log(sg* sqrt(2*3.1416)) + (q2 - 2*mg*q1 + n*mg*mg)/ (sg*sg) + 4  )
}

AIC(49.4347,22.7745,23,,3)

rpsapr = function(sigma,n,t1,t2){
  return(exp( -(9*n/4)*   (1 - ((sqrt((n*t2 - t1^2)/(n^2))) / sigma )^(2/3))^2              ))
}



rpl = function(lambda,x,y){
  return(- (x+y)* log(1+ lambda) + y*log(lambda) + (x+y) *log( 1 + (y/x) )  -  y * log( y/x ))
}


Rpl = function(lambda,x,y){
  return(exp(- (x+y)* log(1+ lambda) + y*log(lambda) + (x+y) *log( 1 + y/x )  -  y * log( y/x )))
}

Rpl1 = function(lambda){
  return(Rpl(lambda,X1,Y1) - 0.1465)
}

Rpl2 = function(lambda){
  return(Rpl(lambda,X2,Y2) - 0.1465)
}

Rpl3 = function(lambda){
  return(Rpl(lambda,X3,Y3) - 0.1465)
}

Rps1 = function(a){
  return(Rps(a,n,t1,t2) - 0.2325)
}
Rps2 = function(a){
  return(Rps(a,n,t1,t2) -  0.1267)
}
Rps3 = function(a){
  return(Rps(a,n,t1,t2) -  0.0289)
}
Rps1b = function(a){
  return(Rps(a,n,q1,q2) - 0.2325)
}
Rps2b = function(a){
  return(Rps(a,n,q1,q2) -  0.1267)
}
Rps3b = function(a){
  return(Rps(a,n,q1,q2) -  0.0289)
}





#-------------------------------------------------


X1 = 49.0366
Y1 = 104.202


X2 = 1762.7062
Y2 = 3999.0244

X3 = 830.5595
Y3 = 1162.3931

plot.function(x = function(t) exp(rpl(t,X1,Y1)),
              from = 0.8,
              to = 4, lwd = 2.5,
              col = "springgreen3",
              main = "Verosimilitud perfil de lambda",
              ylab = "perfil",
              xlab = "Valores")
uno = uniroot(Rpl1, c(1,2))$root
dos = uniroot(Rpl1, c(3,4))$root
segments(x0 = uno, y0 =  0.1465, x1 = dos, y1 =0.1465, col = "steelblue2",
         lwd = 2.5)
abline(v = 1, col = "red", lty=2, lwd=1)
legend("topright", legend="95%",
       lty=1, col = "steelblue2" )
abline(v = Y1 / X1, col = "black", lty=1, lwd=1)
legend("bottomright", legend=expression(paste(hat(lambda))),
       lty=1, col = "black" )


plot.function(x = function(t) exp(rpl(t,X2,Y2)),
              from = 2,
              to = 2.6, lwd = 2.5,
              col = "springgreen3",
              main = "Verosimilitud perfil de lambda2",
              ylab = "perfil",
              xlab = "Valores")
uno = uniroot(Rpl2, c(1,Y2 / X2))$root
dos = uniroot(Rpl2, c(Y2 / X2,3))$root
segments(x0 = uno, y0 =  0.1465, x1 = dos, y1 =0.1465, col = "steelblue2",
         lwd = 2.5)
legend("topleft", legend="95%",
       lty=1, col = "steelblue2" )
abline(v = Y2 / X2, col = "black", lty=1, lwd=1)
legend("bottomright", legend=expression(paste(hat(lambda))),
       lty=1, col = "black" )




plot.function(x = function(t) exp(rpl(t,X3,Y3)),
              from = 1.2,
              to = 1.7, lwd = 2.5,
              col = "springgreen3",
              main = "Verosimilitud perfil de lambda3",
              ylab = "perfil",
              xlab = "Valores")
uno = uniroot(Rpl3, c(1,Y3 / X3))$root
dos = uniroot(Rpl3, c(Y3 / X3,3))$root
segments(x0 = uno, y0 =  0.1465, x1 = dos, y1 =0.1465, col = "steelblue2",
         lwd = 2.5)
abline(v = Y3 / X3, col = "black", lty=1, lwd=1)
legend("bottomright", legend=expression(paste(hat(lambda))),
       lty=1, col = "black" )
legend("topleft", legend="95%",
       lty=1, col = "steelblue2" )



curve(0*x, from = 1, to = 2.4, col = 1)  # Draw Base R plot
curve(exp(rpl(x,X2,Y2)),from = 0, to = 3, col = 1, add = TRUE)
curve(exp(rpl(x,X3,Y3)),from = 0, to = 3, col = 4, add = TRUE)



plot.function(x = function(t) exp(rpl(t,X1,Y1)),
              from = 0.8,
              to = 4, lwd = 2.5,
              col = "black",
              main = "Las 3 verosimilitudes perfiles",
              ylab = "",
              xlab = "")
plot.function(x = function(t) exp(rpl(t,X2,Y2)),
              from = 0.8,
              to = 4, lwd = 2.5,n=10000,
              col = "red",
              main = "Las 3 verosimilitudes perfiles",
              ylab = "",
              xlab = "", add = TRUE)
plot.function(x = function(t) exp(rpl(t,X3,Y3)),
              from = 0.8,
              to = 4, lwd = 2.5,
              col = "blue",,n=10000,
              main = "Las 3 verosimilitudes perfiles",
              ylab = "",
              xlab = "", add = TRUE)
abline(v =Y1 / X1, col = "black", lty=2, lwd=1)
abline(v = Y2 / X2, col = "red", lty=2, lwd=1)
abline(v = Y3 / X3, col = "blue", lty=2, lwd=1)

abline(v = 1, col = "black", lty=1, lwd=1)






#-------------------------Regresión





Datos <- c(0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           0,
           1,
           0,
           1,
           1,
           4,
           0,
           6,
           3,
           8,
           14,
           23,
           3,
           21,
           31,
           34,
           52,
           38,
           34,
           129,
           153,
           156,
           219,
           244,
           17,
           94,
           46,
           27,
           27,
           82,
           43,
           80,
           70,
           31,
           63,
           97,
           90,
           69,
           25,
           146,
           100,
           144,
           92,
           100,
           179,
           253,
           125,
           143,
           166,
           171,
           320,
           269,
           142,
           186,
           248,
           204,
           356,
           299,
           306,
           387,
           450,
           440,
           354,
           237,
           427,
           667,
           528,
           599,
           641,
           702,
           728,
           669,
           790,
           836,
           1167,
           924,
           772,
           808,
           1141,
           994,
           1225,
           1247,
           1038,
           653,
           1683,
           1475,
           1848,
           1737,
           1726,
           1684,
           1464,
           1723,
           3287,
           2658,
           2554,
           2326,
           2610,
           2125,
           2445,
           3166,
           3379,
           3832,
           4328,
           3516,
           2818,
           4103,
           3499,
           3848,
           4996,
           4649,
           4314,
           4545,
           5722,
           6619,
           6253,
           7254,
           6372,
           6167,
           6987,
           8173,
           8781,
           9118,
           10919,
           8826,
           9025,
           10195,
           8863,
           13757,
           12423,
           13579,
           12131,
           11624,
           10559,
           12834,
           13252,
           13454,
           13365,
           13530,
           9356,
           8219,
           13229,
           13183,
           14028,
           12278,
           11301,
           7139,
           7276,
           11431,
           11113,
           11081,
           10168,
           8245,
           5409,
           4483,
           8611,
           8357,
           7336,
           7759,
           6711,
           3762,
           2526,
           2827,
           3970,
           6313,
           4540,
           3714,
           2556,
           2272,
           3940,
           3903,
           3419,
           3729,
           2744,
           1687,
           1576,
           2700,
           2601,
           1857,
           2434,
           2520,
           1997,
           1225,
           2350,
           2435,
           2075,
           1817,
           1643,
           850,
           1086,
           2002,
           2019,
           1972,
           1827,
           1589,
           962,
           777,
           1935,
           2141,
           2067,
           2041,
           1564,
           729,
           1354,
           1918,
           1872,
           1489,
           975,
           1276,
           908,
           908,
           1778,
           1756,
           1760,
           1894,
           1583,
           932,
           1033,
           1925,
           1746,
           1470,
           2559,
           1585,
           893,
           1185,
           1888,
           1781,
           2031,
           1940,
           1672,
           1470,
           1056,
           2067,
           2169,
           1908,
           1845,
           1632,
           896,
           1099,
           1874,
           2068,
           1924,
           1781,
           1379,
           777,
           1249,
           1722,
           1877,
           1772,
           1741,
           1380,
           1255,
           1739,
           2153,
           2352,
           2226,
           2251,
           1853,
           1253,
           1999,
           2905,
           2529,
           3124,
           2662,
           2284,
           2093,
           2508,
           3270,
           3088,
           3390,
           3217,
           2578,
           2316,
           2309,
           4198,
           4427,
           4962,
           4673,
           4141,
           3333,
           4035,
           6750,
           8215,
           8369,
           7930,
           8047,
           5194,
           7598,
           10068,
           9181,
           8778,
           11005,
           9502,
           8842,
           9558,
           14131,
           14391,
           14885,
           11622,
           9559,
           7503,
           9638,
           17817,
           18109,
           16827,
           15093,
           11931,
           12677,
           14497,
           21964,
           21126,
           22113,
           21737,
           17526,
           15137,
           13184,
           18667,
           18615,
           14970,
           14057,
           12341,
           9064,
           9839,
           12787,
           11450,
           11832,
           12345,
           8196,
           4579,
           6078,
           7113,
           7193,
           6178,
           5329,
           4552,
           2563,
           2665,
           4083,
           3774,
           3772,
           3203,
           2450,
           1384,
           1753,
           3178,
           2503,
           2798,
           2396,
           1755,
           1109,
           1217,
           2334,
           2341,
           1923,
           1700,
           1438,
           797,
           1004,
           1873,
           1686,
           1664,
           1456,
           1175,
           569,
           861,
           1456,
           1412,
           1321,
           1234,
           867,
           642,
           997,
           1486,
           1483,
           1232,
           1550,
           1012,
           617,
           939,
           1540,
           1473,
           1471,
           1386,
           1057,
           603,
           513,
           1054,
           1563,
           1525,
           1395,
           971,
           551,
           761,
           1431,
           1302,
           1281,
           782,
           466,
           455,
           440,
           761,
           1374,
           1275,
           1293,
           937,
           659,
           852,
           1609,
           1380,
           1433,
           1333,
           1096,
           748,
           858,
           1578,
           1422,
           1647,
           1393,
           1108,
           854,
           885,
           1258,
           1093,
           1684,
           1642,
           1229,
           902,
           1194,
           2086,
           2162,
           1968,
           2506,
           1789,
           1136,
           1557,
           2776,
           3240,
           3160,
           2767,
           2601,
           1768,
           2369,
           3543,
           3663,
           3352,
           4262,
           2911,
           2397,
           3103,
           4651,
           4449,
           4602,
           4546,
           3778,
           2809,
           3636,
           5817,
           5392,
           5702,
           5483,
           5104,
           3305,
           4234,
           8935,
           9202,
           8068,
           9375,
           7703,
           5582,
           8487,
           13326,
           11838,
           10574,
           13657,
           13235,
           9215,
           11160,
           17599,
           16175,
           18875,
           18065,
           15127,
           12296,
           13428,
           19624,
           21714,
           24417,
           26645,
           16685,
           12589,
           15594,
           21557,
           23048,
           22577,
           21741,
           16401,
           11250,
           12611,
           17595,
           16534,
           16035,
           14790,
           11283,
           7253,
           8982,
           16338,
           14948,
           13802,
           12129,
           9777,
           5701,
           7780,
           17407,
           13708,
           12963,
           12493,
           8783,
           5608,
           9042,
           13343,
           13728,
           13605,
           12384,
           10038,
           6797,
           6580,
           7523,
           14357,
           14000,
           13099,
           10200,
           8031,
           10750,
           14816,
           13754,
           14399,
           13341,
           10813,
           7678,
           10409,
           13331,
           12848,
           12118,
           10234,
           7787,
           5678,
           7127,
           9602,
           9258,
           9255,
           8461,
           5967,
           4143,
           5404,
           7382,
           6307,
           5919,
           5341,
           3985,
           2656,
           3721,
           4695,
           4239,
           3670,
           3306,
           2295,
           1513,
           2210,
           2985,
           2800,
           2275,
           1644,
           973,
           581,
           1375,
           2119,
           1688,
           1645,
           1314,
           814,
           432,
           773,
           1156,
           995,
           930,
           821,
           652,
           202,
           596,
           948,
           953,
           737,
           623,
           417,
           211,
           455,
           595,
           523,
           520,
           452,
           302,
           147,
           333,
           475,
           536,
           478,
           299,
           231,
           107,
           170,
           346,
           321,
           341,
           360,
           206,
           117,
           246,
           307,
           358,
           395,
           308,
           264,
           137,
           275,
           569,
           589,
           794,
           892,
           691,
           314,
           868,
           1275,
           2465,
           2828,
           3220,
           2858,
           2273,
           4373,
           8561,
           11535,
           16055,
           16366,
           11125,
           6381,
           13147,
           19842,
           22388,
           19017,
           17153,
           37875,
           13288,
           23857,
           26389,
           24785,
           20713,
           16080,
           15465,
           8511,
           15423,
           21098,
           21156,
           18847,
           14828,
           5603,
           3778,
           7216,
           9020,
           12978,
           11754
)


plot(Datos, type="l", xlab ="Días", ylab="Casos diarios")


MA <- function(X,n,a){
  m <- c()
  l <- length(X)
  for(i in 1: (l- n+1)){
    
    sum <- 0
    for(j in 1:n){
      sum = sum + X[i+j-1]
    }
    m[i] = sum / n
  }
  
  print(length(seq(n,l,1)) )
  print(length(m) )
  
  if(a==1){
    points(seq(n,l,1) , m, type="l", col = "red")
  }else{
    points(seq((n+1)/2,l - (n-1)/2,1) , m, type="l", col = "red",lwd=2)
  }
  
}
MA(Datos,21,0)







Xvector <- c(30,30,30,40,40,40,40,50,50,50,60,60,60,60,60,70,70,70,70,70)

Xmatriz <- matrix(c(30,30,30,40,40,40,40,50,50,50,60,60,60,60,60,70,70,70,70,70),nrow = 20, ncol = 1)

Yvector <- c(108,110,106,125,120,118,119,132,137,134,148,151,146,147,144,162,156,164,158,159)

Ymatriz = matrix(c(108,110,106,125,120,118,119,132,137,134,148,151,146,147,144,162,156,164,158,159),nrow = 20, ncol = 1)

# Matriz de diseño
K = matrix(c(1,30,1,30,1,30,1,40,1,40,1,40,1,40,1,50,1,50,1,50,1,60,1,60,1,60,1,60,1,60,1,70,1,70,1,70,1,70,1,70),nrow=20,ncol=2,byrow=TRUE)


D = det(t(K) %*% K)

# Tamaño de la muestra
n = 20

# Estimaciones de los parámetros 
bg1 = B1(Xmatriz,Ymatriz,n,D)
bg0 = B0(Xmatriz,Ymatriz,n,bg1)
sg2 = vari(K,Ymatriz,n,bg0,bg1)


# Gráfica de los datos como puntos
plot(Xvector,Yvector,
      pch = 18, ylab = "Presión sistólica",
     xlab = "Edad", main = "Línea de regresión y datos")

# Función de la línea de regresión
regresion = function(t){
  return(bg0 +  bg1*t)
}

plot.function(x = function(t) regresion(t),
              from = 25,
              to = 75, lwd = 2.5,n=2000,
              main = "Línea de regresión y datos",
              col = "red",
              ylab = "Presión sistólica",
              xlab = "Edad", add = TRUE)


# A continuación se grafican las bandas de confianza para la media y los datos

# Esta es la línea inferior de la banda de la media
bandainf1 <- function(x){
  B = matrix(c(bg0,bg1),nrow = 2, ncol=1)
  x0 = matrix(c(1,x),nrow = 1, ncol=2)
  return(bg0 + bg1*x - qt(0.975,18) * sqrt(sg2 * (x0 %*% inv(t(K) %*% K) %*% t(x0))) * sqrt((n-2)/n) )
}


a <- seq(-2.45,80,0.01)

b <- c()

for(i in 1:length(a)){
  b[i] = bandainf1(a[i])
}


points(a,b,type='l',col="blue")

# Esta es la línea superior de la banda de la media

bandasup1 <- function(x){
  B = matrix(c(bg0,bg1),nrow = 2, ncol=1)
  x0 = matrix(c(1,x),nrow = 1, ncol=2)
  return(bg0 + bg1*x + qt(0.975,18) * sqrt(sg2 * (x0 %*% inv(t(K) %*% K) %*% t(x0))) * sqrt((n-2)/n) )
}


a <- seq(-2.45,80,0.01)

b <- c()

for(i in 1:length(a)){
  b[i] = bandasup1(a[i])
}

points(a,b,type='l',col="blue")

# Esta es la línea inferior de la banda de los datos

bandainf2 <- function(x){
  B = matrix(c(bg0,bg1),nrow = 2, ncol=1)
  x0 = matrix(c(1,x),nrow = 1, ncol=2)
  return(bg0 + bg1*x - qt(0.975,18) * sqrt(sg2 * (1+(x0 %*% inv(t(K) %*% K) %*% t(x0)))) * sqrt((n-2)/n) )
}


a <- seq(-2.45,80,0.01)

b <- c()

for(i in 1:length(a)){
  b[i] = bandainf2(a[i])
}

points(a,b,type='l',col="purple")


# Esta es la línea superior de la banda de los datos

bandasup2 <- function(x){
  B = matrix(c(bg0,bg1),nrow = 2, ncol=1)
  x0 = matrix(c(1,x),nrow = 1, ncol=2)
  return(bg0 + bg1*x + qt(0.975,18) * sqrt(sg2 * (1 + (x0 %*% inv(t(K) %*% K) %*% t(x0)))) * sqrt((n-2)/n) )
}


a <- seq(-2.45,80,0.01)

b <- c()

for(i in 1:length(a)){
  b[i] = bandasup2(a[i])
}

points(a,b,type='l',col="purple")




# Residuos de los datos respecto a la línea de regresión
r = Yvector - bg0 - bg1*Xvector


# Gráfica pp de los residuos comparándolos con una distribución normal (0,sg2)
pp_plot(r,0,sqrt(sg2),0.95)



# Logverosimilitud conjunta de b0 y b1
logB <- function(b0,b1){
  B = matrix(c(b0,b1),nrow = 2, ncol=1)
  aux = ( t(Ymatriz - K %*% B)  %*%   (Ymatriz - K %*% B) ) 
  return(-(n* log(aux[1]/n) / 2) - (n/2) )
}

# Contornos de verosimilitud de b0 y b1 con las alturas (0.01,0.05, 0.1)
plotRelative(logB, bg0, bg1, 61, 77, 1.1, 1.5, c(0.01,0.05, 0.1), 500, xlab = "beta 0",ylab = "beta 1")
points( bg0,bg1, type = "p",pch = 8, col = "red")





# Esta función recibe una altura c, una verosimilitud perfil relativa Rp y
# los valores e1,m y e2. La función calcula el intervalo de verosimilitud 
# correspondiente a la altura c en la verosimilitud perfil Rp. Los extremos
# del intervalo se buscarán en los intervalos [e1,m] y [m,e2]
IVc <- function(c,Rp,e1,m,e2){
  
  # Los extremos del intervalo son las raíces de esta función
  f <- function(x){
    return(Rp(x) - c)
  }
  
  uno = uniroot(f, c(e1,m))$root
  dos = uniroot(f, c(m,e2))$root
  return(c(uno,dos))
}





# Estimador de beta 1 restringido a beta 0
EstB1 <- function(b0){
  return((sum(Xvector * Yvector) - b0* sum(Xvector))/sum(Xvector^2))
}

# Logverosimilitud perfil de beta 0
lpb0 <- function(b0){
  res = 0
  for (val in 1: n)
  {
    res = res + (Yvector[val] - b0 - EstB1(b0)*Xvector[val])^2
  }
  return(-(n* log(res/n) / 2) - (n/2))
}

# Verosimilitud perfil relativa de beta 0
Rpb0 <- function(b0){
  return(exp(lpb0(b0) - lpb0(bg0)))
}

# Imprimimos la verosimilitud perfil relativa de beta 0
plot.function(x = function(t) Rpb0(t),
              from = 60,
              to = 78, lwd = 2.5,n=2000,
              main = "Verosimilitud perfil relativa de beta 0",
              ylab = "",
              xlab = "")
abline(v = bg0, col = "red", lty=1, lwd=1)








# Calcularemos las alturas de c correspondientes a los niveles de confianza
# 90, 92, 95 y 99 para beta 0 usando la función IVc y uniroot

# Esta función toma un valor de c, calcula su intervalo de verosimilitud en la 
# relativa perfil de beta 0 y luego regresa la confianza del intervalo 
# segun la cantidad pivotal conocida
confb0 <- function(c){
  a = IVc(c,Rpb0, 0.0001, bg0, 15*bg0)
  return(pt(sqrt((n-2)/ n)*(bg0 - a[1])*sqrt(D/(sg2*sum(Xvector^2)))  ,18) - pt(sqrt((n-2)/ n)*(bg0 - a[2])*sqrt(D/(sg2*sum(Xvector^2)))  ,18))
}

# Estas funciones son tales que sus raíces estiman las alturas c para las confianzas
# 90%, 92%, 95% y 99% respectivamente .
Estc90 <- function(c){
  return(confb0(c) - 0.90)
}

Estc92 <- function(c){
  return(confb0(c) - 0.92)
}

Estc95 <- function(c){
  return(confb0(c) - 0.95)
}

Estc99 <- function(c){
  return(confb0(c) - 0.99)
}


# Se estiman las alturas para las distintas confianzas.
c95 = uniroot(Estc95, c(0.0001,1))$root
c99 = uniroot(Estc99, c(0.0001,1))$root
c90 = uniroot(Estc90, c(0.0001,1))$root
c92 = uniroot(Estc92, c(0.0001,1))$root


# Dadas las alturas, calculamos numéricamente los intervalos asociadas a éstas
# en la verosimilitud perfil de beta 0. Luego graficamos cada intervalo
f90 <- function(x){
  return(Rpb0(x) - c90)
}
uno = uniroot(f90, c(bg0 / 10,bg0))$root
dos = uniroot(f90, c(bg0,10*bg0))$root
segments(x0 = uno, y0 = c90, x1 = dos, y1 = c90, col = "steelblue2",
         lwd = 2.5)

f92 <- function(x){
  return(Rpb0(x) - c92)
}
uno = uniroot(f92, c(bg0 / 10,bg0))$root
dos = uniroot(f92, c(bg0,10*bg0))$root
segments(x0 = uno, y0 = c92, x1 = dos, y1 = c92, col = "deeppink3",
         lwd = 2.5)


f95 <- function(x){
  return(Rpb0(x) - c95)
}
uno = uniroot(f95, c(bg0 / 10,bg0))$root
dos = uniroot(f95, c(bg0,10*bg0))$root
segments(x0 = uno, y0 = c95, x1 = dos, y1 = c95, col = "chartreuse4",
         lwd = 2.5)

f99 <- function(x){
  return(Rpb0(x) - c99)
}
uno = uniroot(f99, c(bg0 / 10,bg0))$root
dos = uniroot(f99, c(bg0,10*bg0))$root
segments(x0 = uno, y0 = c99, x1 = dos, y1 = c99, col = "purple",
         lwd = 2.5)
legend("topright", legend=c("90%", "92%","95%","99%"),
       lty=1, col = c("steelblue2", "deeppink3" ,"chartreuse4" , "purple") )






# Estimador de beta 0 restringido a beta 1
EstB0 <- function(b1){
  return((sum(Yvector) - b1* sum(Xvector))/n)
}

# Logverosimilitud perfil de beta 1
lpb1 <- function(b1){
  res = 0
  for (val in 1: 20)
  {
    res = res + (Yvector[val] - EstB0(b1) - b1*Xvector[val])^2
  }
  return(-(n* log(res/n) / 2) - (n/2))
}

# Verosimilitud perfil relativa de beta 1
Rpb1 <- function(b1){
  return(exp(lpb1(b1) - lpb1(bg1)))
}

# Graficamos la verosimilitud perfil relativa de beta 1
plot.function(x = function(t) Rpb1(t),
              from = 1.1,
              to = 1.5, lwd = 2.5,n=2000,
              main = "Verosimilitud perfil relativa de beta 1",
              ylab = "",
              xlab = "")
abline(v = bg1, col = "red", lty=1, lwd=1)


# Calcularemos las alturas de c correspondientes a los niveles de confianza
# 90, 92, 95 y 99 para beta 1 usando la función IVc y uniroot
confb1 <- function(c){
  a = IVc(c,Rpb1, bg1/15, bg1, 15*bg1)
  return(pt(sqrt((n-2)/ n)*(bg1 - a[1])*sqrt(D/(n*sg2))  ,18) - pt(sqrt((n-2)/ n)*(bg1 - a[2])*sqrt(D/(sg2*n))  ,18))
}

# Estas funciones son tales que sus raíces estiman las alturas c para las confianzas
# 90%, 92%, 95% y 99% respectivamente .
Estc90 <- function(c){
  return(confb1(c) - 0.90)
}

Estc92 <- function(c){
  return(confb1(c) - 0.92)
}

Estc95 <- function(c){
  return(confb1(c) - 0.95)
}

Estc99 <- function(c){
  return(confb1(c) - 0.99)
}


# Se estiman las alturas para las distintas confianzas.
c95 = uniroot(Estc95, c(0.0001,1))$root
c99 = uniroot(Estc99, c(0.0001,1))$root
c90 = uniroot(Estc90, c(0.0001,1))$root
c92 = uniroot(Estc92, c(0.0001,1))$root


# Dadas las alturas, calculamos numéricamente los intervalos asociadas a éstas
# en la verosimilitud perfil de beta 1. Luego graficamos cada intervalo
f90 <- function(x){
  return(Rpb1(x) - c90)
}
uno = uniroot(f90, c(bg1 / 15,bg1))$root
dos = uniroot(f90, c(bg1,300*bg1))$root
segments(x0 = uno, y0 = c90, x1 = dos, y1 = c90, col = "steelblue2",
         lwd = 2.5)

f92 <- function(x){
  return(Rpb1(x) - c92)
}
uno = uniroot(f92, c(bg1 / 10,bg1))$root
dos = uniroot(f92, c(bg1,10*bg1))$root
segments(x0 = uno, y0 = c92, x1 = dos, y1 = c92, col = "deeppink3",
         lwd = 2.5)


f95 <- function(x){
  return(Rpb1(x) - c95)
}
uno = uniroot(f95, c(bg1 / 10,bg1))$root
dos = uniroot(f95, c(bg1,10*bg1))$root
segments(x0 = uno, y0 = c95, x1 = dos, y1 = c95, col = "chartreuse4",
         lwd = 2.5)

f99 <- function(x){
  return(Rpb1(x) - c99)
}
uno = uniroot(f99, c(bg1 / 10,bg1))$root
dos = uniroot(f99, c(bg1,10*bg1))$root
segments(x0 = uno, y0 = c99, x1 = dos, y1 = c99, col = "purple",
         lwd = 2.5)
legend("topright", legend=c("90%", "92%","95%","99%"),
       lty=1, col = c("steelblue2", "deeppink3" ,"chartreuse4" , "purple") )











# Verosimilitud perfil relativa de sigma^2, recibe el valor a evaluar, el vector Y,
# la matriz de diseño K y las estimaciones b0, b1 y sg correspondientes a beta 0, beta 1
# y sigma^2
Rps2 <- function(sigma2,Y,K,b0,b1,sg){
  B = matrix(c(b0,b1),nrow = 2, ncol=1)
  aux = ( t(Y - K %*% B)  %*%   (Y - K %*% B) ) 
  return(exp((-(n* log(sigma2))/2 - (aux[1]) / (2*sigma2) ) -   (-(n* log(sg))/2 - (aux[1]) / (2*sg) )        )       )
} 

# Esta función es la verosimilitu perfil relativa de sigma^2 pero ya con los parámetros
# del ejemplo. Se hace ésta para el posterior cálculo de los intervalos de confianza con
# la cantidad pivotal.
Rpsv <- function(t){
  return(Rps2(t,Ymatriz,K,bg0,bg1,sg2))
}

# Graficamos la verosimilitud perfil relativa de beta 1
plot.function(x = function(t) Rps2(t,Ymatriz,K,bg0,bg1,sg2),
              from = 2,
              to = 20, lwd = 2.5,
              main = "Verosimilitud perfil relativa de sigma cuadrada",
              ylab = "",
              xlab = "")
abline(v = sg2, col = "red", lty=1, lwd=1)





# Calcularemos las alturas de c correspondientes a los niveles de confianza
# 90, 92, 95 y 99 para sigma^2 usando la función IVc y uniroot
confsigma2 <- function(c){
  a = IVc(c,Rpsv,0.0001,sg2,15*sg2)
  return(pchisq(n*sg2 / a[1] ,18) - pchisq(n*sg2 / a[2] ,18))
}

# Estas funciones son tales que sus raíces estiman las alturas c para las confianzas
# 90%, 92%, 95% y 99% respectivamente .
Estc90 <- function(c){
  return(confsigma2(c) - 0.90)
}

Estc92 <- function(c){
  return(confsigma2(c) - 0.92)
}

Estc95 <- function(c){
  return(confsigma2(c) - 0.95)
}

Estc99 <- function(c){
  return(confsigma2(c) - 0.99)
}


# Se estiman las alturas para las distintas confianzas.
c95 = uniroot(Estc95, c(0.0001,1))$root
c99 = uniroot(Estc99, c(0.0001,1))$root
c90 = uniroot(Estc90, c(0.0001,1))$root
c92 = uniroot(Estc92, c(0.0001,1))$root


# Dadas las alturas, calculamos numéricamente los intervalos asociadas a éstas
# en la verosimilitud perfil de sigma^2. Luego graficamos cada intervalo.
f90 <- function(x){
  return(Rpsv(x) - c90)
}
uno = uniroot(f90, c(0.0001,sg2))$root
dos = uniroot(f90, c(sg2,15*sg2))$root
segments(x0 = uno, y0 = c90, x1 = dos, y1 = c90, col = "steelblue2",
         lwd = 2.5)

f92 <- function(x){
  return(Rpsv(x) - c92)
}
uno = uniroot(f92, c(0.0001,sg2))$root
dos = uniroot(f92, c(sg2,15*sg2))$root
segments(x0 = uno, y0 = c92, x1 = dos, y1 = c92, col = "deeppink3",
         lwd = 2.5)


f95 <- function(x){
  return(Rpsv(x) - c95)
}
uno = uniroot(f95, c(0.0001,sg2))$root
dos = uniroot(f95, c(sg2,15*sg2))$root
segments(x0 = uno, y0 = c95, x1 = dos, y1 = c95, col = "chartreuse4",
         lwd = 2.5)

f99 <- function(x){
  return(Rpsv(x) - c99)
}
uno = uniroot(f99, c(0.0001,sg2))$root
dos = uniroot(f99, c(sg2,15*sg2))$root
segments(x0 = uno, y0 = c99, x1 = dos, y1 = c99, col = "purple",
         lwd = 2.5)
legend("topright", legend=c("90%", "92%","95%","99%"),
       lty=1, col = c("steelblue2", "deeppink3" ,"chartreuse4" , "purple") )






