function read_emis,file

t=create_struct('freq',           0.0, $
                'sst',            0.0, $
                'wsp',            0.0, $
                'emis_rss', fltarr(2), $
                'emis_fm4', fltarr(2), $
                'emis_fm5', fltarr(2), $
                'emis_fm6', fltarr(2))

a=replicate(t,16,16,6)

openr,1,file
readu,1,a
close,1

return,a
end
