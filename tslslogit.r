## Copyright 2018 Zach Flynn

## Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

## The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

tslslogit <- function (y, r, w=r, tol=1e-8, noconstant=FALSE)
{
  ## Model is:
  ## E[y*w] = E[(1+exp(-1*r'beta))^(-1)*w]

  if (ncol(r)!=ncol(w))
    stop("Number of covariates must equal number of instruments")

  if (!noconstant)
  {
    r <- cbind(1,r)
    w <- cbind(1,w)
  }
  
  wr <- t(r) %*% w
  wy <- t(w) %*% y

  qr.wr <- qr(wr)

  gamma.hat <-  backsolve(qr.R(qr.wr), t(qr.Q(qr.wr)) %*% wy)

  beta <- rep(0,times=ncol(r))

  err <- 1
  while (err > tol)
  {
    wybeta <- t(w) %*% (1/(1+exp(-1*r %*% beta)))
    gamma.beta <- backsolve(qr.R(qr.wr), t(qr.Q(qr.wr)) %*% wybeta)
    err <- max(abs(gamma.beta - gamma.hat))
    beta <- beta + gamma.hat - gamma.beta
  }
   
  return (beta)
}

