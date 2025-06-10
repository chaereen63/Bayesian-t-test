#### asht package source code (central Behrens-Fisher distribution) ####
# 원본 pbefi 함수 (paste.txt에서 가져온 것)
pbefi<-function(q,n1,n2,R=NULL,s1=NULL,s2=NULL,epsilon=10^(-8)){
  if (is.null(R) & (is.numeric(s1) & is.numeric(s2))){
    R<-atan((s1/sqrt(n1))/(s2/sqrt(n2)))
  } else if (is.null(s1) & is.null(s2)){
    # do nothing R<-R
  } else stop("supply one of R or (s1 and s2)")
  
  sinR<-sin(R)
  cosR<-cos(R)
  
  # if sinR=0 then cosR=1 and D=T2
  # and if cosR=0 then sinR=1 and D=T1
  if (abs(sinR)<epsilon){
    out<-pt(q,n2-1)
  } else if (abs(cosR)<epsilon){
    out<-pt(q,n1-1)
  } else if (cosR>0 & cosR<1 & sinR>0 & sinR<1){
    ## usual case       
    ifunc<-function(u,Y){
      pt((Y-u*sinR)/cosR,n2-1)*dt(u,n1-1)
    }
    nq<-length(q)
    out<-rep(NA,nq)
    for (i in 1:nq){
      out[i]<-integrate(ifunc,-Inf,Inf,Y=q[i])$value
    }
  }
  
  out
}

# PDF 계산을 위한 함수 (수치 미분)
dbefi <- function(x, n1, n2, s1=NULL, s2=NULL, R=NULL, h=0.001) {
  (pbefi(x+h, n1, n2, s1=s1, s2=s2, R=R) - pbefi(x-h, n1, n2, s1=s1, s2=s2, R=R)) / (2*h)
}

# 시각화 함수
plot_befi_central <- function(n1, n2, s1, s2, xlim=c(-4, 4)) {
  
  # R 매개변수 계산
  R <- atan((s1/sqrt(n1))/(s2/sqrt(n2)))
  
  # x 값들
  x <- seq(xlim[1], xlim[2], length.out=200)
  
  # PDF와 CDF 계산
  pdf_vals <- sapply(x, function(t) dbefi(t, n1, n2, s1=s1, s2=s2))
  cdf_vals <- pbefi(x, n1, n2, s1=s1, s2=s2)
  
  # 비교를 위한 t-분포
  t1_pdf <- dt(x, n1-1)
  t2_pdf <- dt(x, n2-1)
  t1_cdf <- pt(x, n1-1)
  t2_cdf <- pt(x, n2-1)
  
  # 플롯 설정 - 2x1 레이아웃으로 더 크게
  par(mfrow=c(1,2), mar=c(4,4,3,2), mgp=c(2.5,1,0))
  
  # PDF 플롯
  plot(x, pdf_vals, type='l', lwd=3, col='blue', 
       main=paste("Behrens-Fisher Distribution\nn1=", n1, ", n2=", n2, ", s1=", s1, ", s2=", s2),
       xlab='x', ylab='Density', ylim=c(0, max(pdf_vals, t1_pdf, t2_pdf)*1.05))
  lines(x, t1_pdf, col='red', lty=2, lwd=2)
  lines(x, t2_pdf, col='darkgreen', lty=3, lwd=2)
  
  legend('topright', 
         legend=c('Behrens-Fisher', paste('t(', n1-1, ')', sep=''), paste('t(', n2-1, ')', sep='')),
         col=c('blue', 'red', 'darkgreen'), lty=c(1,2,3), lwd=c(3,2,2),
         cex=0.8)
  
  # CDF 플롯
  plot(x, cdf_vals, type='l', lwd=3, col='blue',
       main="Cumulative Distribution Function",
       xlab='x', ylab='P(X ≤ x)', ylim=c(0,1))
  lines(x, t1_cdf, col='red', lty=2, lwd=2)
  lines(x, t2_cdf, col='darkgreen', lty=3, lwd=2)
  
  legend('bottomright', 
         legend=c('Behrens-Fisher', paste('t(', n1-1, ')', sep=''), paste('t(', n2-1, ')', sep='')),
         col=c('blue', 'red', 'darkgreen'), lty=c(1,2,3), lwd=c(3,2,2),
         cex=0.8)
  
  # 매개변수 정보 출력
  cat("\n=== 매개변수 정보 ===\n")
  cat("n1 =", n1, ", n2 =", n2, "\n")
  cat("s1 =", s1, ", s2 =", s2, "\n")
  cat("R =", round(R, 4), "(", round(R*180/pi, 1), "°)\n")
  cat("s1/√n1 =", round(s1/sqrt(n1), 4), "\n")
  cat("s2/√n2 =", round(s2/sqrt(n2), 4), "\n")
  
  # 분위수 비교 테이블
  quantiles <- c(0.025, 0.05, 0.95, 0.975)
  
  # BEFI 분위수 계산
  befi_quant <- sapply(quantiles, function(p) {
    uniroot(function(x) pbefi(x, n1, n2, s1=s1, s2=s2) - p, 
            interval=c(-10, 10))$root
  })
  
  t1_quant <- qt(quantiles, n1-1)
  t2_quant <- qt(quantiles, n2-1)
  
  cat("\n=== 분위수 비교 ===\n")
  cat(sprintf("%-8s %-12s %-10s %-10s\n", "분위수", "BEFI", paste("t(",n1-1,")",sep=""), paste("t(",n2-1,")",sep="")))
  cat(sprintf("%-8s %-12s %-10s %-10s\n", "------", "----------", "-------", "-------"))
  
  for(i in 1:length(quantiles)) {
    cat(sprintf("%-8s %-12.3f %-10.3f %-10.3f\n", 
                paste(quantiles[i]*100, "%", sep=""),
                befi_quant[i], t1_quant[i], t2_quant[i]))
  }
}

#### noncentral Behrens-Fisher distribution ####
# 비중심 Behrens-Fisher CDF 함수
pbefi_noncentral <- function(q, n1, n2, delta=0, R=NULL, s1=NULL, s2=NULL, epsilon=10^(-8)) { # delta = mu1-mu2
  if (is.null(R) & (is.numeric(s1) & is.numeric(s2))){
    R<-atan((s1/sqrt(n1))/(s2/sqrt(n2)))
  } else if (is.null(s1) & is.null(s2)){
    # do nothing R<-R
  } else stop("supply one of R or (s1 and s2)")
  
  sinR<-sin(R)
  cosR<-cos(R)
  
  if (abs(sinR)<epsilon){
    # 비중심 t-분포
    ncp <- delta / sqrt(s2^2/n2)
    out<-pt(q, n2-1, ncp=ncp)
  } else if (abs(cosR)<epsilon){
    # 비중심 t-분포
    ncp <- delta / sqrt(s1^2/n1)
    out<-pt(q, n1-1, ncp=ncp)
  } else if (cosR>0 & cosR<1 & sinR>0 & sinR<1){
    # 일반적인 경우: 적분식에 delta 추가
    ifunc<-function(u,Y){
      pt((Y - delta - u*sinR)/cosR, n2-1) * dt(u, n1-1)
    }
    nq<-length(q)
    out<-rep(NA,nq)
    for (i in 1:nq){
      out[i]<-integrate(ifunc,-Inf,Inf,Y=q[i])$value
    }
  }
  out
}

# 비중심 BEFI PDF (수치 미분)
dbefi_noncentral <- function(x, n1, n2, delta=0, s1=NULL, s2=NULL, R=NULL, h=0.001) {
  (pbefi_noncentral(x+h, n1, n2, delta, s1=s1, s2=s2, R=R) - 
     pbefi_noncentral(x-h, n1, n2, delta, s1=s1, s2=s2, R=R)) / (2*h)
}

# 비중심 Behrens-Fisher 분포 시각화 함수
plot_befi_noncentral <- function(n1, n2, s1, s2, delta=0, xlim=c(-4, 4)) {
  
  # R 매개변수 계산
  R <- atan((s1/sqrt(n1))/(s2/sqrt(n2)))
  
  # x 값들
  x <- seq(xlim[1], xlim[2], length.out=200)
  
  # PDF와 CDF 계산 (비중심)
  pdf_vals <- sapply(x, function(t) dbefi_noncentral(t, n1, n2, delta=delta, s1=s1, s2=s2))
  cdf_vals <- pbefi_noncentral(x, n1, n2, delta=delta, s1=s1, s2=s2)
  
  # 비교를 위한 중심 분포
  pdf_central <- sapply(x, function(t) dbefi(t, n1, n2, s1=s1, s2=s2))
  cdf_central <- pbefi(x, n1, n2, s1=s1, s2=s2)
  
  # 비교를 위한 t-분포
  t1_pdf <- dt(x, n1-1)
  t2_pdf <- dt(x, n2-1)
  t1_cdf <- pt(x, n1-1)
  t2_cdf <- pt(x, n2-1)
  
  # 플롯 설정 - 여백을 더 타이트하게 조정
  par(mfrow=c(1,2), 
      mar=c(4, 4, 3, 1),    # 오른쪽 여백을 1로 줄임
      mgp=c(2.5, 1, 0),
      oma=c(0, 0, 0, 0))    # 외부 여백 제거
  
  # PDF 플롯
  plot(x, pdf_vals, type='l', lwd=3, col='blue', 
       main=paste("Noncentral Behrens-Fisher Distribution\nn1=", n1, ", n2=", n2, ", s1=", s1, ", s2=", s2, ", δ=", delta),
       xlab='x', ylab='Density', 
       ylim=c(0, max(pdf_vals, pdf_central, t1_pdf, t2_pdf)*1.05))
  lines(x, pdf_central, col='blue', lty=2, lwd=2)
  lines(x, t1_pdf, col='red', lty=3, lwd=2)
  lines(x, t2_pdf, col='darkgreen', lty=4, lwd=2)
  
  # 범례 위치와 크기 조정
  legend('topleft', # 범례 위치
         legend=c(paste('BEFI (δ=', delta, ')', sep=''), 'BEFI (δ=0)', 
                  paste('t(', n1-1, ')', sep=''), paste('t(', n2-1, ')', sep='')),
         col=c('blue', 'blue', 'red', 'darkgreen'), 
         lty=c(1,2,3,4), lwd=c(3,2,2,2),
         cex=0.65,              # 범례 텍스트 크기 더 작게
         bg='white',            # 배경색 추가
         box.lwd=0.5,           # 테두리 선 얇게
         inset=c(0.02, 0.02))   # 가장자리에서 약간 떨어뜨림
  
  # CDF 플롯
  plot(x, cdf_vals, type='l', lwd=3, col='blue',
       main="Cumulative Distribution Function",
       xlab='x', ylab='P(X ≤ x)', ylim=c(0,1))
  lines(x, cdf_central, col='blue', lty=2, lwd=2)
  lines(x, t1_cdf, col='red', lty=3, lwd=2)
  lines(x, t2_cdf, col='darkgreen', lty=4, lwd=2)
  
  # 범례 위치와 크기 조정
  legend('bottomright', 
         legend=c(paste('BEFI (δ=', delta, ')', sep=''), 'BEFI (δ=0)', 
                  paste('t(', n1-1, ')', sep=''), paste('t(', n2-1, ')', sep='')),
         col=c('blue', 'blue', 'red', 'darkgreen'), 
         lty=c(1,2,3,4), lwd=c(3,2,2,2),
         cex=0.65,              # 범례 텍스트 크기 더 작게
         bg='white',            # 배경색 추가
         box.lwd=0.5,           # 테두리 선 얇게
         inset=c(0.02, 0.02))   # 가장자리에서 약간 떨어뜨림
  
  # 매개변수 정보 출력
  cat("\n=== 매개변수 정보 ===\n")
  cat("n1 =", n1, ", n2 =", n2, "\n")
  cat("s1 =", s1, ", s2 =", s2, "\n")
  cat("delta (δ) =", delta, "\n")
  cat("R =", round(R, 4), "(", round(R*180/pi, 1), "°)\n")
  cat("s1/√n1 =", round(s1/sqrt(n1), 4), "\n")
  cat("s2/√n2 =", round(s2/sqrt(n2), 4), "\n")
  
  # 분위수 비교 테이블
  quantiles <- c(0.025, 0.05, 0.95, 0.975)
  
  # 비중심 BEFI 분위수
  befi_nc_quant <- sapply(quantiles, function(p) {
    tryCatch({
      uniroot(function(x) pbefi_noncentral(x, n1, n2, delta=delta, s1=s1, s2=s2) - p, 
              interval=c(-20, 20))$root
    }, error = function(e) NA)
  })
  
  # 중심 BEFI 분위수
  befi_c_quant <- sapply(quantiles, function(p) {
    uniroot(function(x) pbefi(x, n1, n2, s1=s1, s2=s2) - p, 
            interval=c(-10, 10))$root
  })
  
  t1_quant <- qt(quantiles, n1-1)
  t2_quant <- qt(quantiles, n2-1)
  
  cat("\n=== 분위수 비교 ===\n")
  cat(sprintf("%-8s %-10s %-10s %-10s %-10s\n", "분위수", "BEFI(δ≠0)", "BEFI(δ=0)", paste("t(",n1-1,")",sep=""), paste("t(",n2-1,")",sep="")))

  for(i in 1:length(quantiles)) {
    cat(sprintf("%-8s %-10.3f %-10.3f %-10.3f %-10.3f\n", 
                paste(quantiles[i]*100, "%", sep=""),
                befi_nc_quant[i], befi_c_quant[i], t1_quant[i], t2_quant[i]))
  }
  
}


# plot_befi_noncentral(n1=30, n2=30, s1=2, s2=2, delta=1.5)
