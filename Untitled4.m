function out = fft_paolo(X,N,PI)    
% called by : RAYGANS, PHO   
% Computes Fast Fourier Transform
%      COMPLEX*16 X(1),U,W,T
% COMPLEX*8 X(1),U,W,T

M=log(N)/log(2.)+.1;
NV2=N/2;
NM1=N-1;
J=1;

for I=1:NM1
    if (I >= J) GO TO 10
      T=X(J)
      X(J)=X(I)
      X(I)=T
10    K=NV2
20    IF(K .GE. J) GO TO 30
      J=J-K
      K=K/2
      GO TO 20
30    J=J+K
40    CONTINUE
%
      DO 70 L=1,M
      LE=2**L
      LE1=LE/2
      U=(1.,0.)
      PF=PI/FLOAT(LE1)
      W=CMPLX(COS(PF),-SIN(PF))
      DO 60 J=1,LE1
      DO 50 I=J,N,LE
      IP=I+LE1
      T=X(IP)*U
      X(IP)=X(I)-T
      X(I)=X(I)+T
50    CONTINUE
      U=U*W
60    CONTINUE
70    CONTINUE
%
%              NORMALIZATION
      X(1)=X(1)/N
      DO 80 I=2,N
80    X(I)=X(I)*2/N
      RETURN
      END