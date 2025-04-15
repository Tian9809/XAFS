clear all
close all

disp(' ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp('%         CONTINUOUS CAUCHY WAVELET TRANSFORM     % ')
disp('%                 OF EXAFS SIGNAL                 % ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp(' ')
disp('code freely downloaded from http://www.univ-mlv.fr/~farges/waw ')
disp('(c) 2000, Univ. Marne la Vallee, France ')
disp(' ')
disp('Our reference paper for this code is : ')
disp('  Munoz M., Argoul P. et Farges F. (2003) ')
disp('  Continuous Cauchy wavelet transform analyses of EXAFS spectra: a qualitative approach.  ')
disp('  American Mineralogist volume 88, pp. 694-700. ')
disp(' ')


% version history:
% 1999 Hans-Argoul : core wavelet algorithm
% 1999-2002 Argoul-Munoz : EXAFS adapation
% 2002 Farges : graphical and user interface
% 2003 Munoz : CPU optimizations and graphical updates
% 2003 Farges-Munoz : various fixes and web version


%--------------------------------------------------------------------
% USER PREFERENCES : CHECK THESE BEFORE RUNNING SCRIPT

% mode : 1 for computation+plot     2 for plot only (an already computed .m file)
mode = 1 ;

% Cauchy order (recommended for standard analysis: 200)
n = 200 ;

% minimum R-space distance ? (recommended: 0.2; 0 is forbidden)
ri = 0.2 ;

% maximum R-space distance ? (recommended: 6.0)
rf = 6 ;

% number of R-space intervals (recommended: 200)
na = 200 ;

%palette for the 2D CCWT : 0 (no 2D display)   1 (color)   2 (grays)
palette1 = 1 ;

%palette for the 3D CCWT : 0 (no 3D display)   1 (color)   2 (grays)   3 (grays/satin look)
palette2 = 3 ;

%viewing angle for the 3D : 0 (defaut view)   1  OR  2  OR 3   (various other angles: view(-45,65), view(-125,55) OR view(44,44))
vieww = 3 ;
%--------------------------------------------------------------------



% code starts here

% case compute-and-plot

if mode == 1 
%%%%%%%%%%%%%%%%%%%%%%
%  Input ascii file  %
%%%%%%%%%%%%%%%%%%%%%%

[fichier, path] = uigetfile('*', 'Enter normalized EXAFS spectrum: ');
fid = fopen(fichier,'r');
[A,count] = fscanf(fid,'%g %g',[2 inf]);
A = A';
kold = A(:,1);
xold = A(:,2);
status = fclose(fid);

% removing extension quand besoin est
if findstr(fichier,'.')
    fichier = fichier(1:findstr(fichier,'.')-1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAFS data interpolation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('EXAFS data interpolation...')
nt = 256;
kfin = kold(length(kold));
pask = (kfin-kold(1))/nt;
knew = [kold(1):pask:kfin]';
xnew = interp1q(kold,xold,knew);
tab1 = [knew xnew];
clear tab1 chi kfin kold nt tab1 xold
knew = knew'; xnew = xnew';


  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Wavelet transform analysis       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FT parameters
z=2^3;                                  
nk = length(knew);
ZF = z*nk;			            		
npt = ZF/2;			            		
freq = 1/pask*(0:npt-1)/ZF;	    		
omega = 2*pi*freq;		

% TF calculation
tff = fft(xnew(1,:),ZF); 
TF = tff(1:ZF/2)./max(knew);

% Scale parameter
pasr = (rf-ri)/na;							  
for i = 1:(na+1) 			                  
  r(i) = ri + (i-1)*pasr;					  
  a(i) = n/2/r(i);							 
end

% Characteristic values of the Cauchy wavelet
derpha0 = n+1;                                
s = sum(log(1:n));
maxi = exp(log(2*pi) + n*log(n) - n - s);     


% Main calculation
% Cauchy wavelet calculation
disp('     ')
disp('Cauchy wavelet calculation...')
disp('     ')
for i = 1:(na+1)                  	
  int = a(i).*omega;
  for j = 1:npt
    if int(j) == 0
      filtre(i,j) = 0.;
    else
      filtre(i,j) = exp(log(2*pi) - s + n*log(int(j)) - int(j));
    end
  end
end
clear int

% Wavelet transform calculation
disp('Wavelet transform calculation...')
disp('     ')
for i = 1:(na+1)       				
  noyau = conj(filtre(i,:)).*tff(1:npt);
  res = ifft(noyau,ZF);		  
  to(1:nk,i) = res(1:nk).';
end
clear filtre noyau


% save TO
 fichier2 = strcat(fichier,'.mat');
 wt_data = abs(to);
 save (fichier2,'knew','xnew','freq','TF','ri','rf','to','a','n', 'wt_data') 
 
else

  
% mode 2 : entree du nom du fichier .mat
  [fichier2, path] = uigetfile('*.mat', 'Enter *.mat file : ');
  load (fichier2,'knew','xnew','freq','TF','ri','rf','to','a','n')
  
  % enlever extension
  fichier = fichier2(1:findstr(fichier2,'.mat')-1)
  
end





% PLOTS
if palette1 ~= 0
     
figure('Position',[1, 100, 950, 650])
%clf
set(gcf,'Color','white')

% 1er s/plot : TF
axes('position', [0.07, 0.35, 0.21, 0.62])
plot(abs(TF),freq.*pi,'k','linewidth',1);
set(gca,'XDir','reverse','YDir','normal');
set(gca,'FontName','Helvetica')
set(gca,'FontSize',14)
xlabel('FT-modulus','FontName','Helvetica','FontSize',18);
ylabel( strcat('R + {\Delta}R (',setstr( hex2dec('c5')),')'  ),'FontName','Helvetica','FontSize',18);
axis([0 max(abs(TF))*1.05 ri rf]);

% 2eme s/plot : WT
axes('position', [0.3, 0.35, 0.67, 0.62])
pcolor(knew,n/2./a,abs(to)') 
axis off
shading interp

if palette1 == 1
   colormap(jet)
 else
   colormap(1-pink)
end

set(gca,'FontName','Helvetica') 
set(gca,'FontSize',14)

% 3eme s/plot : EXAFS
axes('position', [0.3, 0.11, 0.67, 0.22]);
plot(knew,xnew,'k','linewidth',1);
zz = max(abs(min(xnew)),max(xnew));
set(gca,'FontName','Helvetica')
set(gca,'FontSize',14)
axis([min(knew) max(knew) min(xnew)*1.05 max(xnew)*1.05]);
xlabel( strcat('{\itk (}', setstr( hex2dec('c5') ), '^-^1)' ), 'FontName','Helvetica','FontSize',18); 
ylabel('\chi ({\itk})', 'FontName','Helvetica','FontSize',18);

% sauvegarde fichier format JPEG 100%
  fichier3 = strcat(fichier,'_TO2D');
  print( gcf, '-djpeg100', fichier3 )

end



% 3D graph of CCWT modulus

 if palette2 ~= 0
  figure(2)
  set(gcf,'Color','white')
  
  if palette2 == 1
     surf(knew,n/2./a,abs(to)');
     colormap(jet)
 end

  if palette2 == 2
     surf(knew,n/2./a,abs(to)');
     colormap(1-pink)
 end
  
 if palette2 == 3
     colormap(1-gray)
     surfl(knew,n/2./a,abs(to)');
 end
  
  shading interp
  grid off
  axis([min(knew) max(knew) min(n/2./a) max(n/2./a) ]);
   
  set(gca,'FontName','Helvetica') 
  set(gca,'FontSize',18)
  xlabel(strcat('{\itk (}', setstr( hex2dec('c5') ), '^-^1)' ) ); 
  ylabel(strcat('R + {\Delta}R (',setstr( hex2dec('c5')),')'  ) );
  zlabel('CCWT modulus')
  
  if vieww == 1
    view(-45,65)
  end
  if vieww == 2
    view(-125,55)
  end
  if vieww == 3
     view(44,44)
  end
  
  rotate3d on
   
   
  % affiche la legende des couleurs
  colorbar
  
 
  % sauvegarde fichier format JPEG 100%
  fichier3 = strcat(fichier,'_TO3D');
  print( gcf, '-djpeg100', fichier3 )
 
end


% the end
disp('End of computation, have a good day.')
disp('      ')