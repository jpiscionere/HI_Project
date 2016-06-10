compile_opt idl2, strictarrsubs
base =  '/data0/LasDamas/Consuelo/'
filebase = '/data0/LasDamas/Consuelo/SEED/Consuelo_SEED_z0p054_fof_b0p2.fdpp.halos'
start_seed = 1
end_seed =  50
seed_offset = 4000
Nbins =  10
Npartmin =  500.0d
Npartmax =  2500.0d
savfile = 'allnpart.sav'
if n_elements(allnpart) eq 0 then begin
   if (findfile(savfile))[0] ne '' then begin
      restore, savfile, /verbose
   endif else begin
      allnpart =  ptrarr(end_seed-start_seed+1)
   endelse
endif

allmf =  dblarr(end_seed-start_seed+1, Nbins)
for seed = start_seed, end_seed do begin
   file = str_replace( filebase, 'SEED', strn(seed+seed_offset), /global)
   print, file, format = '("Reading file `",A0,"`..",$)'
   if ptr_valid(allnpart[seed-start_seed]) eq 0 then begin
      rdfloat, file, junk, npart, /double, /silent
      ;;    readcol, file,  npart, format = 'X,L', /silent ;;; skip column 1
      allnpart[seed-start_seed] =  ptr_new(npart, /no_copy)
   endif
   npart =  *allnpart[seed-start_seed]
   hist =  double(histogram( alog10(npart), min = alog10(Npartmin), max = alog10(Npartmax), Nbins = Nbins))   
   allmf[seed-start_seed, 0:Nbins-1] =  hist
   print, '. done'
endfor

;;; compute the median mass function
median_mf =  median(allmf, dimension = 1)
binsize =  alog10(Npartmax/Npartmin)/Nbins
mass =  10.0d^(dindgen(Nbins)*binsize + alog10(Npartmin))

size = 800
position = [0.2, 0.4, 0.9, 0.9]
residual_position = [0.2, 0.2, 0.9, 0.4]
xrange = [npartmin, npartmax]
if !d.name ne 'PS' then $
   window, xsize = size, ysize = size

plot, mass, median_mf, /ylog, /xlog, xrange = xrange, position = position, $
      xtickformat = '(A1)'
wanted_seeds =  [42, 1, 36, 34, 38]
wanted_seeds =  [49, 15, 31, 26, 38]
;wanted_seeds =  [42, 36]
;wanted_seeds = lindgen(50)+1
nseeds =  n_elements(wanted_seeds)
deviations =  dblarr(nseeds, Nbins)
sum_deviations =  dblarr(nseeds)
alltot_signs = fltarr(nseeds)
for iseed = 0, nseeds-1 do begin
   seed =  wanted_seeds[iseed]
   this_mf = allmf[seed-start_seed, 0:Nbins-1]
   if nseeds lt 6 then $
      oplot, mass, this_mf, color = collist[iseed mod n_elements(collist)], psym = -cgsymcat((iseed+4))
   deviation =  reform((this_mf - median_mf)/sqrt(median_mf))
   deviations[iseed, 0:Nbins-1] = deviation
   sum_deviations[iseed] =  total(deviation^2, /preserve)
   signs =  sign(deviation)
   alltot_signs[iseed]  =  total(signs)
endfor
if nseeds lt 6 then $
   legend, string(wanted_seeds, format = '(I0)')+'='+strcompress(string(sum_deviations, format = '(F0.2)'), /remove_all), textcolor = collist[0:nseeds-1], /top, /right, box = 0

yrange = [-2.0, 2.0]
plot, [0], /nodata, /noerase, xrange = xrange, yrange = yrange, $
      position = residual_position, ytickinterval = 1.0, yminor = 1, /xlog
plots, xrange, [0.0, 0.0], line = 2, thick = 3

for iseed = 0, nseeds-1 do begin
   deviation =  deviations[iseed, 0:Nbins-1]
   if nseeds lt 6 then begin
      oplot, mass, deviation, color = collist[iseed mod n_elements(collist)], psym = -cgsymcat((iseed+4))
      forprint, mass, deviation, /text
   endif
endfor
xx =  bsort(sum_deviations)
forprint, wanted_seeds[xx], sum_deviations[xx], alltot_signs[xx], format = '(I3," ",F12.2," ",I3)', /text
stop
;forprint, wanted_seeds[xx], sum_deviations[xx], format = '(I3," ",F12.2)', textout = 'Consuelo_z0.0_mf_ranks(nearM_star).txt',/nocomment



;;; Compute the D-statistic from K-S test assuming the truth is the
;;; median mass function.
d_stat =  dblarr(end_seed-start_seed+1)
prob_seed =  dblarr(end_seed-start_seed+1)
for seed = start_seed, end_seed do begin
   this_mf = allmf[seed-start_seed, 0:Nbins-1]
   kstwo, median_mf, this_mf, D, prob
   d_stat[seed-start_seed] = D
   prob_seed[seed-start_seed] =  prob
endfor
allseeds =  lindgen(end_seed-start_seed+1)+start_seed
xx =  bsort(d_stat)
d_stat =  d_stat[xx]
prob_seed =  prob_seed[xx]
allseeds =  allseeds[xx]
forprint, allseeds, d_stat,prob_seed /text

stop



end
