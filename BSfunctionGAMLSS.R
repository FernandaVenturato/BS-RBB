#----------------------------------------------------------------------------------------
# Utilizando a distribuicao BS em GAMLSS
#----------------------------------------------------------------------------------------
BS <- function(mu.link = "log", sigma.link = "log")
{
  mstats <- checklink("mu.link", "BS", substitute(mu.link),
                      c("inverse", "log", "identity", "own"))
  dstats <- checklink("sigma.link", "BS", substitute(sigma.link),
                      c("inverse", "log", "identity", "own"))
  structure(
    list(family = c("BS", "Birnbaum-Saunders"),
         parameters = list(mu=TRUE,sigma=TRUE),
         nopar = 2,
         type = "Continuous",
         mu.link = as.character(substitute(mu.link)),
         sigma.link = as.character(substitute(sigma.link)),
         mu.linkfun = mstats$linkfun,
         sigma.linkfun = dstats$linkfun,
         mu.linkinv = mstats$linkinv,
         sigma.linkinv = dstats$linkinv,
         mu.dr = mstats$mu.eta,
         sigma.dr = dstats$mu.eta,
         dldm = function(y,mu,sigma) 1/mu^3*(y/sigma+sigma/y-2)-1/mu,
         d2ldm2 = function(mu) -2/mu^2,
         dldd = function(y,mu,sigma) 1/(y+sigma)-1/(2*sigma)-1/(2*mu^2)*(1/y-y/sigma^2),
         d2ldd2 = function(sigma) -(1/(y+sigma)-1/(2*sigma)-1/(2*mu^2)*(1/y-y/sigma^2))^2,
         d2ldmdd = function(y) rep(0,length(y)),
         G.dev.incr = function(y,mu,sigma,...) -2*dBS(y,mu,sigma,log=TRUE),
         rqres = expression(rqres(pfun="pBS", type="Continuous",y=y,mu=mu,sigma=sigma)),
         mu.initial = expression({ mu <- (y+mean(y))/2 }),
         sigma.initial = expression({sigma <- rep(median(y),length(y))}),
         mu.valid = function(mu) all(mu > 0) ,
         sigma.valid = function(sigma) all(sigma > 0),
         y.valid = function(y) all(y > 0)
    ),
    class = c("gamlss.family","family"))
}

#----------------------------------------------------------------------------------------
# fdp da distribuicao BS
dBS<-function(x, mu=1, sigma=1, log=FALSE)
{
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  if (any(x < 0)) stop(paste("x must be positive", "\n", ""))
  fy <- 1/(2*mu*sqrt(2*pi*sigma))*x^(-3/2)*(x+sigma)*exp(-1/(2*mu^2)*(x/sigma+sigma/x-2))
  fy
}
#----------------------------------------------------------------------------------------
# fda da distribuicao BS
pBS <- function(q, mu=1, sigma=1, lower.tail = TRUE, log.p = FALSE)
{
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  if (any(q < 0)) stop(paste("y must be positive", "\n", ""))
  cumulative <- function(x, a, b)
  {
    integration <- integrate(dBS, lower = 0, upper = x, mu = a, sigma = b)$value
    return(integration)
  }
  cdf <- mapply(cumulative, q, a = mu, b = sigma)
  if(lower.tail==FALSE) cdf <- 1-cdf
  if(log.p==TRUE) cdf <- log(cdf)
  cdf
}
#----------------------------------------------------------------------------------------
# funcao quantilica da distribuicao BS
qBS <- function(p, mu=1, sigma=1, lower.tail = TRUE, log.p = FALSE)
{
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  if (log.p==TRUE) p <- log(p)
  if (lower.tail==FALSE) p <- 1-p
  if (any(p < 0)|any(p > 1)) stop(paste("p must be between 0 and 1", "\n", ""))
  q <- sigma*(mu*qnorm(p)/2+sqrt((mu*qnorm(p)/2)^2)+1)^2
  return(q)
}
#----------------------------------------------------------------------------------------
# gerador de numeros aleatorios da distribuicao BS
rBS <- function(n, mu=1, sigma=1)
{
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  if (any(n <= 0)) stop(paste("n must be a positive integer", "\n", ""))
  n <- ceiling(n)
  p <- runif(n)
  r <- qBS(p,mu=mu,sigma=sigma)
  r
}
#----------------------------------------------------------------------------------------