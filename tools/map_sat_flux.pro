;-----------------------------------------------------------------------------
; 
;                                 map_sat_flux.pro
;
; Program reads data from *_satellite_*.flux and plots flux on a
; geographic map.
;
; Created on 23 September 2024 by Mei-Ching Fok (Code 673/GSFC)
;-----------------------------------------------------------------------------

  close,/all

; get colors for color table palette
  loadct,39                               ;  rainbow with white
  tvlct, red, green, blue, /get
  !p.background=255
  !p.color=0
  !p.charthick=3  

; Setup color range 
  colormax=250
  colormin=10 
  ncolor=colormax-colormin

; Read flux size and setup arrays 
  header=' '
  filehead=' '
  read,'enter filehead of satellite flux (i.e., 2017a147_storie_e) => ',filehead
  openr,1,filehead+'.flux' 
  readf,1,npos,je
  glon=fltarr(npos)   &  glat=glon    &   Lshell=glon
  gride=fltarr(je)     &  flux1=gride
  fluxa=fltarr(npos,je)
  Elab=strarr(je+1)  &  Elabi=Elab
  readf,1,header
  readf,1,gride      ; energy in keV
  for i=1,2 do readf,1,header

; Read satellite location and flux
  for n=0,npos-1 do begin
      readf,1,iyr,iday,ihr,imin,sec,Ls
      readf,1,xgeo,ygeo,zgeo
      readf,1,flux1
      rgeo=sqrt(xgeo*xgeo+ygeo*ygeo+zgeo*zgeo)
      glat(n)=asin(zgeo/rgeo)*180./!pi
      glon(n)=atan(ygeo,xgeo)*180./!pi
      Lshell(n)=Ls
      fluxa(n,*)=flux1(*)
  endfor
  close,1
  ntime=npos

; make 2D arrays of glon2D, glat2D, and L2D
  jpos=94      ; no. of minutes in one orbit
  iorb=20
  glon2D=fltarr(iorb,jpos)  &  glat2D=glon2D  &  L2D=glon2D
  for i=0,iorb-1 do begin
      n1=i*jpos
      n2=n1+jpos-1
      glon2D(i,0:jpos-1)=glon(n1:n2)
      glat2D(i,0:jpos-1)=glat(n1:n2)
      L2D(i,0:jpos-1)=Lshell(n1:n2)
  endfor

; Setup energy boundaries (ebound) and find flux in energy range
  ebound=fltarr(je+1)
  fluxk=fltarr(ntime)
  for k=1,je-1 do ebound(k)=sqrt(gride(k-1)*gride(k))
  ebound(0)=gride(0)*gride(0)/ebound(1)
  ebound(je)=gride(je-1)*gride(je-1)/ebound(je-1)
  for k=0,je do begin
      Elab(k)=string(ebound(k),'(f6.1)')
      Elabi(k)=string(round(ebound(k)),'(i4.4)')
  endfor
  Ehead=' '
  print,'plot energy(keV) : '
  for k=0,je-1 do print,k,'  '+Elab(k),' - ',Elab(k+1)
  read,' enter lower and upper energy bins => ',ie1,ie2
  Ehead=Elabi(ie1)+'-'+Elabi(ie2+1)+'keV'
  for n=0,ntime-1 do begin
      fluxk(n)=0.
      for k=ie1,ie2 do fluxk(n)=fluxk(n)+fluxa(n,k)*(ebound(k+1)-ebound(k))
  endfor

; Setup flux range    
  yon=' '
  fmax=alog10(max(fluxk))
  print,'max log flux is => ',fmax  
  read,'Do you want to change it ? (y)yes, (n)no => ',yon
  if (yon eq 'y') then read,'enter new max log flux => ',fmax
  fmin=fmax-3.     
  maxV=string(fmax,'(f4.1)')
  minV=string(fmin,'(f4.2)')
  cof=ncolor/(fmax-fmin)

; Setup for plotting flux and color bars
  xi=0.07
  xf=0.9
  y0=0.1 
  y2=0.8
  dy=(y2-y0)/ncolor
  yave=0.5*(y0+y2)
  bar_lab='log flux (#/cm2/s/sr)'

; Plot satellite flux
  set_plot,'ps'
  device,filename=filehead+'_'+Ehead+'.ps',/inches,yoffset=3.,$
             xoffset=0.3,xsize=8.,ysize=5.36,/color,bits_per_pixel=24
  MAP_SET,/CYLINDRICAL,/GRID,/CONTINENTS,/LABEL,title=filehead+'  '+Ehead, $
          pos=[xi,y0,xf,y2],limit=[-90,-180,90,180],/noerase
  fmin1=10.^fmin
  for n=0,ntime-1 do begin
      logflux=fmin
      if (fluxk(n) gt fmin1) then logflux=alog10(fluxk(n))
      if (logflux gt fmax) then logflux=fmax
      fcolor=colormin+(logflux-fmin)*cof    
      oplot,[glon(n),glon(n)],[glat(n),glat(n)],color=fcolor,psym=1
  endfor
  contour,L2D,glon2D,glat2D,levels=[2,3],thick=2,/overplot

; plot color bar
  x1=xf+0.04
  x2=x1+0.02
  for i=0,ncolor do begin
      y1=y0+i*dy
      polyfill,[x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],color=fix(i),/normal
  endfor
  xyouts,x2+0.01,y0,minV,/normal
  xyouts,x2+0.01,y1-0.01,maxV,/normal
  xyouts,x1-0.01,yave,bar_lab,orientation=90,alignment=0.5,/normal
     
  device,/close_file
  end
