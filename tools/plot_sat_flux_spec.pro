;-------------------------------------------------------------------------------
;
;                           plot_sat_flux_spec.pro
;
; IDL procedure reads *.flux file and plot the energy-Time spectrogram or 
; line plot
;
; Created on 14 April 2015 by Mei-Ching Fok, Code 673, NASA GSFC.
;-------------------------------------------------------------------------------

close,/all

; Setup plot  
x_wsize=600
y_wsize=300

; get colors for color table palette
ibar=1
;read,'color bar? (1)rainbow, (2)HENA, (3)Akebono => ',ibar
if (ibar eq 1) then begin
   loadct,39       ;  rainbow with white
   tvlct, red, green, blue, /get
endif else begin
   red=fltarr(256)    &    green=fltarr(256)    &    blue=fltarr(256)
   if (ibar eq 2) then begin
      openr,1,'hena.tbl'
      readf,1,red,green,blue
   endif else begin
      openr,1,'akebono.tbl'
      for i=0,255 do begin
          readf,1,ic,red1,green1,blue1
          red(i)=red1
          green(i)=green1
          blue(i)=blue1
      endfor
      blue(0:1)=[0.,1.]
   endelse
   close,1
   tvlct, red, green, blue
endelse
!p.background=255
!p.charthick=4
!p.thick=3

; read file name
fhead=' '
read,'enter the file head (e.g., 2013a076_rbspA_e) => ',fhead
openr,2,fhead+'.flux'
read,'hour label? (1) decimal hour, (2) hhmm => ',topt
x_tit='hour'
if (topt eq 2) then x_tit='hhmm'
ig=0
readf,2,ntime,je         ; je=number of energy grid, ig=number of pa grid
hour=fltarr(ntime)  &  Lshell=hour  
flux=fltarr(ntime,je+6)  &  flux1=fltarr(je)   &  energy=flux1  & Penergy=flux1
fluxPe=flux   &   fluxPa=flux      ; perpendicular and parallel flux
ratio=flux
header=' '
readf,2,header 
readf,2,energy       ; energy in keV
if (ig gt 0) then begin
   pitA=fltarr(ig)
   readf,2,header
   readf,2,pitA
endif
readf,2,header 
bx=0. &  by=0.  &  bz=0.
ratio(*,*)=0.
for n=0,ntime-1 do begin
    readf,2,year,month,day,hr,minu,sec,L1,density,Efield
    readf,2,flux1
    for k=0,je-1 do flux(n,k)=flux1(k)*energy(k)
    readf,2,flux1
    for k=0,je-1 do fluxPe(n,k)=flux1(k)*energy(k)     ; Perpendicular flux
    readf,2,flux1
    for k=0,je-1 do fluxPa(n,k)=flux1(k)*energy(k)     ; Parallel flux
    for k=0,je-1 do begin
        if (fluxPe(n,k) gt 1.e3) then ratio(n,k)=fluxPa(n,k)/fluxPe(n,k)    
    endfor
    flux(n,je)=density/1.e6       ; density in cm-3
    flux(n,je+5)=Efield
    if (n eq 0) then day0=day
    if (topt eq 1) then hour(n)=(day-day0)*24.+hr+minu/60.+sec/3600.
    if (topt eq 2) then hour(n)=hr*100+minu+sec/60.
    Lshell(n)=L1
    for m=0,ig-1 do readf,2,flux1
endfor 
close,2
date=string(fix(month),'(i2.2)')+string(fix(day0),'(i2.2)')

; Remove singular points in plasmasphere density
for n=1,ntime-1 do begin
    fof=flux(n,je)/flux(n-1,je)
    if (fof lt 1.e-3) then flux(n,je)=flux(n-1,je)
endfor

; setup species
nsp=strlen(fhead)-1
spe=strmid(fhead,nsp,1)                ; species
if (spe eq 'h') then species='H+'
if (spe eq 'o') then species='O+'
if (spe eq 'e' or spe eq '1') then species='e-'

new_plot:
; Choose between energy-time spectrogram or line plot
iplt=1
;read,'plot? (1)E-Time-Spectro (2)Line-Plot (3)plasma-den (4)B (5)E => ',iplt
if (iplt eq 1) then read, $
   ' plot? (1)omni (2)perpendicular (3) parallel-flux (4) jpar/jper => ',iopt
ihead=' '
if (iopt eq 1) then ihead=' j_omni'
if (iopt eq 2) then ihead=' j_perp'
if (iopt eq 3) then ihead=' j_para'
if (iopt eq 4) then ihead=' jpar/jper'
yon=' '
icolor=0
ie1=0
ie2=fix(je-1)
if (iplt eq 2) then begin
   print, 'energy bins from ie1, ie2 => ',ie1,ie2
   read,'Do you want them? y(yes), n(no) => ',yon
   if (yon eq 'y') then read,'enter new ie1,ie2 => ',ie1,ie2
   read,'dark color for (1)low-energy, (2)high-energy => ',icolor
endif
if (iplt eq 4) then begin
   ie1=je+1
   ie2=je+4
   icolor=2
endif
if (iplt eq 5) then begin
   ie1=je+5
   ie2=je+5
   icolor=2
endif
iunit=-1
ieu=0
if (iplt le 2) then begin  
  if (iopt lt 4) then read, $
     'unit (0)flux/(cm2ssr) (1)flux/(keVcm2ssr) (2)flux/(MeVcm2ssr) => ',iunit
  read,'energy in? (0) MeV, (1) keV, (2) eV => ',ieu
endif
read,'plot in?  0 = linear scale, or 1 = log scale => ',ilog
plog=0
if (iopt le 3 and ilog eq 1) then plog=1     ; label color bar with log value

Eunit='keV'
Penergy=energy
if (ieu eq 0) then begin
   Eunit='MeV'
   Penergy=energy/1000.    ; energy in MeV
endif
if (ieu eq 2) then begin
   Eunit='eV'
   Penergy=energy*1000.    ; energy in eV
endif
if (iopt eq 1) then fluxs=flux
if (iopt eq 2) then fluxs=fluxPe
if (iopt eq 3) then fluxs=fluxPa
if (iopt eq 4) then fluxs=ratio 
if (iunit eq 2) then fluxs(*,*)=fluxs(*,*)*1000.
if (iopt eq 4) then flx_lab=ihead
if (iunit eq 0) then flx_lab='eflux(keV/cm2/sr/s)'
if (iunit eq 1) then flx_lab='eflux(keV/cm2/sr/s/keV)'
if (iunit eq 2) then flx_lab='eflux(keV/cm2/sr/s/MeV)'
if (iplt eq 1) then y_title='energy ('+Eunit+')'
if (iplt eq 2) then y_title=flx_lab        
if (iplt eq 3) then y_title='density (cm-3)'
if (iplt eq 4) then y_title='Bgsm (nT)'
if (iplt eq 5) then y_title='E field (mV/m)'
if (plog eq 1) then flx_lab='log '+flx_lab

; setup plot ranges
if (iplt le 2) then fmax=max(fluxs(*,0:je-1))    
if (iplt eq 3) then fmax=max(fluxs(*,je))    
if (iplt eq 4) then fmax=max(fluxs(*,je+4))    
if (iplt eq 5) then fmax=max(fluxs(*,je+5))    
print,' fmax = ',fmax
read,' Do you want to change it? (y/n) => ',yon
if (yon eq 'y') then read,' enter new fmax => ',fmax
fmin=0.   
if (ilog eq 1 and iplt eq 1) then fmin=fmax/1.e5
if (iplt eq 4) then fmin=-fmax
print,' fmin = ',fmin
read,' Do you want to change it? (y/n) => ',yon
if (yon eq 'y') then read,' enter new fmin => ',fmin
if(min(fluxs) lt fmin) then fluxs(where(fluxs lt fmin))=fmin
pmin=fmin
pmax=fmax
if (plog eq 1) then pmin=alog10(fmin)
if (plog eq 1) then pmax=alog10(fmax)
Emin=Penergy(0)
Emax=Penergy(je-1)
if (iplt eq 1) then begin
 ; print,'Emin, Emax = ',Emin,Emax
 yon='n'
 ; read,' Do you want to change it? (y/n) => ',yon
   if (yon eq 'y') then read,' enter new Emin, Emax => ',Emin,Emax
endif

set_plot,'ps'
file_name=fhead+'_flux.ps'
if (iplt eq 3) then file_name=fhead+'_PlDen.ps'
if (iplt eq 4) then file_name=fhead+'_B.ps'
if (iplt eq 5) then file_name=fhead+'_E.ps'
device,filename=file_name,/inches,/color,yoffset=3.0, $
       xoffset=0.3,xsize=8.,ysize=8.*y_wsize/x_wsize

; setup levels
nlevel=1
if (iplt eq 1) then nlevel=59
if (iplt eq 2 or iplt eq 4) then nlevel=ie2-ie1+1
lvl=fltarr(nlevel)    &    colr=intarr(nlevel)
dlvl=(fmax-fmin)/(nlevel-1)
if (ilog eq 1) then dlvl=alog10(fmax/fmin)/(nlevel-1)
colrmax=254           
colrmin=1
ncolor=colrmax-colrmin+1
dcolr=(float(colrmax)-float(colrmin))/(nlevel-1)
for i=0,nlevel-1 do begin
    if (ilog eq 0) then lvl(i)=fmin+i*dlvl
    if (ilog eq 1) then begin
       lvl1=alog10(fmin)+i*dlvl
       lvl(i)=10.^lvl1
    endif
    colr(i)=round(float(colrmin)+i*dcolr)
endfor
if (iplt eq 1) then begin
   ilow=2
 ; read,'color flux beyond lower limit (1) black, (2) white => ',ilow
   if (ilow eq 2) then colr(0)=255    
endif

; Plot plasma flux
x0=0.16
xf=x0+0.7
y0=0.16 
yf=0.81
ym=0.5*(y0+yf)
thr=24
hr_max=ceil(max(hour)/thr)*thr
hr_min=fix(min(hour)/thr)*thr
x_tick=(hr_max-hr_min)*2/thr
hr_max=max(hour)
hr_min=min(hour)
x_tick=ceil((hr_max-hr_min)*60.)
print,'hr_min,hr_max,x_tick => ',hr_min,hr_max,x_tick
read,' Do you want to change it? (y/n) => ',yon
if (yon eq 'y') then read,'new hr_min,hr_max,x_tick => ',hr_min,hr_max,x_tick
x_minor=(hr_max-hr_min)/x_tick
if (iplt eq 1) then begin
   contour,fluxs(*,0:je-1),hour,Penergy,xtitle=x_tit, $
        ytitle=y_title,xrange=[hr_min,hr_max],yrange=[Emin,Emax],xstyle=1, $
        ystyle=1,xticks=x_tick,xminor=x_minor,title=fhead+ihead,charsize=1.0, $
        pos=[x0,y0,xf,yf],color=0,levels=lvl,c_colors=colr,ylog=1,/fill 
 ; oplot,hour,fluxs(*,je),color=230,thick=3
endif else begin
   k1=je
   if (icolor eq 1) then k1=ie1
   if (icolor eq 2) then k1=ie2  
   plot,hour,fluxs(*,k1),xtitle=x_tit,ytitle=y_title,yrange=[fmin,fmax], $
        xrange=[hr_min,hr_max],xstyle=1,ystyle=1,xticks=x_tick,title=fhead,$
        charsize=2.0,pos=[x0,y0,xf,yf],color=colr(0)
   if (iplt eq 4) then oplot,[hr_min,hr_max],[0.,0.],linestyle=1
   for k=ie1,ie2 do begin
       if (icolor eq 1) then k1=k-ie1
       if (icolor eq 2) then k1=ie2-k
       if (iplt eq 2 or iplt eq 4) then oplot,hour,fluxs(*,k),color=colr(k1)
   endfor
endelse
xyouts,x0,yf+0.02,date,size=1.4,/normal

; draw color bar and label
plab=' '
x1=xf+0.04
x2=x1+0.03
if (iplt eq 1) then dy=(yf-y0)/(colrmax-colrmin+1)
if (iplt eq 2 or iplt eq 4) then dy=(yf-y0)/nlevel
if (iplt eq 1) then begin
   for i=colrmin,colrmax do begin
       y1=y0+i*dy
       y2=y1+dy*1.03
       polyfill,[x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],color=i,/normal  
   endfor
   xyouts,x2,y0-0.02,string(pmin,'(f4.1)'),size=1.4,color=0,/normal
   xyouts,x2,yf,string(pmax,'(f4.1)'),size=1.4,color=0,/normal
   xyouts,x2+0.03,ym,flx_lab,alignment=0.5,orientation=90.,/normal 
endif
if (iplt eq 2 or iplt eq 4) then begin
   y1=y0+0.5*dy
   for k=ie2,ie1,-1 do begin
      if (icolor eq 1) then k1=k-ie1
      if (icolor eq 2) then k1=ie2-k
      if (iplt eq 2) then begin
         if (iunit eq 0) then plab='>'+string(Penergy(k),'(f5.0)')+' '+Eunit 
         if (iunit ne 0) then plab=string(Penergy(k),'(f5.1)')+' '+Eunit 
      endif else begin
         if (k eq je+1) then plab='Bx'
         if (k eq je+2) then plab='By'
         if (k eq je+3) then plab='Bz'
         if (k eq je+4) then plab='B'
      endelse
      xyouts,xf+0.02,y1,plab,color=colr(k1),size=1.5,/normal
      y1=y1+dy
   endfor
endif   
    
device,/close_file
read,'Do you want to continue?  (y)yes, (n)no => ',yon
if (yon eq 'y') then goto,new_plot

end
