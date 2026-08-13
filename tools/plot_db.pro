 close,/all

 loadct,39       ;  rainbow with white
 tvlct, red, green, blue, /get
 !p.color=0
 !p.background=255
 !p.charsize=2
 !p.charthick=3  
 !p.thick=3  
 DEVICE, DECOMPOSED=0

 header=' '
 fhead=' ' 
 read,' enter the file head of .db (i.e., 2013d151) => ', fhead
 read,' plot: 0=*.dbe, 1=*.dbi, 2=*.db, 3=*.RBC => ',isp
 read,'Enter number of time => ',np    
 itime=1     ; 1 = hour, 2 = day
 ftail=' '
 if (isp eq 0) then ftail='dbe'   
 if (isp eq 1) then ftail='dbi'   
 if (isp eq 2) then ftail='db'   
 if (isp eq 3) then ftail='RBC'  
 openr,1,fhead+'.'+ftail

 sps='(i)'
 if (isp eq 0) then sps='(e)'
 if (isp eq 2) then sps='(i+e)'
 if (isp eq 3) then sps='(RBC)'
 hour=fltarr(np)   &   Dst=hour  &  DstRC=hour  &  totE=hour  &  totEB=hour
 day=hour
 readf,1,header
 for m=0,np-1 do begin
     if (isp lt 3) then readf,1,hour1,Kp,AE,AL,Dst1,DstRC1,rcsum,totE1,totEB1
     if (isp eq 3) then readf,1,hour1,Kp,AE,AL,Dst1,DstRC1,totE1,totEB1
     hour(m)=hour1
     day(m)=hour1/24.
     Dst(m)=Dst1
     DstRC(m)=DstRC1
     totE(m)=totE1
     totEB(m)=totEB1
 endfor
 close,1

 Dstmin=-150.  
 Dstmax=40.
 Emin=0.
 Emax=3.0e31
 if (isp eq 2) then Emax=4.e31
 if (isp eq 3) then Emax=1.e28

 tmin=0.
 daymax=ceil(max(day))
 hourmax=daymax*24
 if (itime eq 1) then begin
    time=hour
    time_lab='hour'
    tmax=hourmax
 endif else begin
    time=day 
    time_lab='day'
    tmax=daymax
 endelse

 set_plot,'ps'
 device,filename=fhead+'_'+ftail+'.ps'

 iplt=1
 iplt_lab=' '
 if (iplt eq 1) then iplt_lab='-B '
 if (iplt eq 1) then totEP=totEB
 if (iplt eq 2) then totEP=totE
 if (isp eq 3) then totEP=totE
 pos1=[0.1,0.12,0.84,0.9]
 plot,time,Dst,title=fhead,xtitle=time_lab,ytitle='Dst (nT)', $
      xrange=[tmin,tmax],xstyle=1,pos=pos1,ystyle=1, $
      yrange=[Dstmin,Dstmax],linestyle=1      ; dotted line
 oplot,time,DstRC,linestyle=2                 ; dashed line
 axis,yaxis=1,yrange=[Dstmin,Dstmax],ystyle=1,color=255
 plot,time,totEP,xrange=[tmin,tmax],xstyle=1,pos=pos1, $
      yrange=[Emax,Emin],ystyle=4,/noerase
 axis,yaxis=1,ytitle='Ring Current Energy (keV)'
 xyouts,0.5,0.85,'..... Dst',/normal
 xyouts,0.5,0.80,'-- Dst*',/normal
 xyouts,0.5,0.75,'___ Erc'+iplt_lab+sps,/normal

 device,/close_file
 end
