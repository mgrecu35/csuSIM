pro plt_emis,vmin=vmin,vmax=vmax,vpol=vpol,hpol=hpol,ps=ps,chs=chs

if (n_elements(vmin) EQ 0) then vmin=-10
if (n_elements(vmax) EQ 0) then vmax=10
if (n_elements(chs)  EQ 0) then chs=0.7

t=create_struct('freq',           0.0, $
                'sst',            0.0, $
                'wsp',            0.0, $
                'emis_rss', fltarr(2), $
                'emis_fm4', fltarr(2), $
                'emis_fm5', fltarr(2), $
                'emis_fm6', fltarr(2))

a=replicate(t,16,16,6)

if (keyword_set(vpol)) then outfile='cmpr_emis-vpol.ps' else outfile='cmpr_emis-hpol.ps'
openr,1,'cmpr_emis.dat'
readu,1,a
close,1

for i=0,5 do begin
  freq = a[0,0,i].freq
  ttl1='(FASTEM4 - FASTEM6)'
  ttl2='(FASTEM5 - FASTEM6)'
  ttl3='(RSS - FASTEM6)'
  dtv_fm4 = a[*,*,i].sst * (a[*,*,i].emis_fm4[0] - a[*,*,i].emis_fm6[0])
  dth_fm4 = a[*,*,i].sst * (a[*,*,i].emis_fm4[1] - a[*,*,i].emis_fm6[1])
  dtv_fm5 = a[*,*,i].sst * (a[*,*,i].emis_fm5[0] - a[*,*,i].emis_fm6[0])
  dth_fm5 = a[*,*,i].sst * (a[*,*,i].emis_fm5[1] - a[*,*,i].emis_fm6[1])
  dtv_rss = a[*,*,i].sst * (a[*,*,i].emis_rss[0] - a[*,*,i].emis_fm6[0])
  dth_rss = a[*,*,i].sst * (a[*,*,i].emis_rss[1] - a[*,*,i].emis_fm6[1])

  if (keyword_set(vpol)) then begin
    ttl=string(freq,format='("Freq = ",f6.2," (V-Pol)")')
    pimage,dtv_fm4,position=ppos(6,3,i+0),/multi,col_file='dcolor.tbl',proj='xy',xlabel='wsp (m/sec)', $
                   ylabel='sst (K)',region=[275,305,0,30],vmin=vmin,vmax=vmax,xws=2000,yws=1200,title=ttl, $
                   subtitle=ttl1,ps_file=outfile,ps_land=ps,ls=chs,ts=chs,units='Ts Diff (K)'
    pimage,dtv_fm5,position=ppos(6,3,i+6),/multi,col_file='dcolor.tbl',proj='xy',xlabel='wsp (m/sec)', $
                   ylabel='sst (K)',region=[275,305,0,30],vmin=vmin,vmax=vmax,title=ttl,subtitle=ttl2, $
                   ps_land=ps,ls=chs,ts=chs,units='Ts Diff (K)'
    pimage,dtv_rss,position=ppos(6,3,i+12),/multi,col_file='dcolor.tbl',proj='xy',xlabel='wsp (m/sec)', $
                   ylabel='sst (K)',region=[275,305,0,30],vmin=vmin,vmax=vmax,title=ttl,subtitle=ttl3, $
                   ps_land=ps,ls=chs,ts=chs,units='Ts Diff (K)'
  endif else begin
    ttl=string(freq,format='("Freq = ",f6.2," (H-Pol)")')
    pimage,dth_fm4,position=ppos(6,3,i+0),/multi,col_file='dcolor.tbl',proj='xy',xlabel='wsp (m/sec)', $
                   ylabel='sst (K)',region=[275,305,0,30],vmin=vmin,vmax=vmax,xws=2000,yws=1200,title=ttl, $
                   subtitle=ttl1,ps_file=outfile,ps_land=ps,ls=chs,ts=chs,units='Ts Diff (K)'
    pimage,dth_fm5,position=ppos(6,3,i+6),/multi,col_file='dcolor.tbl',proj='xy',xlabel='wsp (m/sec)', $
                   ylabel='sst (K)',region=[275,305,0,30],vmin=vmin,vmax=vmax,title=ttl,subtitle=ttl2, $
                   ps_land=ps,ls=chs,ts=chs,units='Ts Diff (K)'
    pimage,dth_rss,position=ppos(6,3,i+12),/multi,col_file='dcolor.tbl',proj='xy',xlabel='wsp (m/sec)', $
                   ylabel='sst (K)',region=[275,305,0,30],vmin=vmin,vmax=vmax,title=ttl,subtitle=ttl3, $
                   ps_land=ps,ls=chs,ts=chs,units='Ts Diff (K)'
  endelse
endfor
!p.multi[0]=0
if (keyword_set(ps)) then begin
  device,/close_file
  set_plot,'x'
endif
end
