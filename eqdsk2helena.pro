; Convert a G-EQDSK format file into input for HELENA
;
;

PRO eqdsk2helena, input, psin=psin, mfm=mfm, nrmap=nrmap, npmap=npmap, new=new
  
  safe_colors, /first

  IF NOT KEYWORD_SET(psin) THEN psin=0.995
  IF NOT KEYWORD_SET(mfm) THEN mfm = 64
  IF NOT KEYWORD_SET(nrmap) THEN nrmap = 101
  IF NOT KEYWORD_SET(npmap) THEN npmap = 129

  IF NOT is_pow2(mfm) THEN BEGIN
    PRINT, "ERROR: mfm must be a power of 2"
    RETURN
  ENDIF
  
  IF NOT is_pow2(npmap-1) THEN BEGIN
    PRINT, "ERROR: npmap must be 2^n + 1"
    RETURN
  ENDIF

  ; Read the grid input
  grid = read_neqdsk(input)
   
  psi_boundary = grid.sibdry
  psi_axis = grid.simagx
  
  ; Contour psi to get the plasma boundary
  
  contour_lines, grid.psi, findgen(grid.nx), findgen(grid.ny), $
    levels=[psi_axis + (psi_boundary - psi_axis)*psin], $
    path_info=info, path_xy=xy
  
  IF N_ELEMENTS(info) GT 1 THEN BEGIN
    ; Find the surface closest to the magnetic axis

    ; Need the index of the magnetic axis
    ri = INTERPOL(findgen(grid.nx), grid.r[*,0], grid.rmagx)
    
    zi = INTERPOL(findgen(grid.ny), grid.z[0,*], grid.zmagx)

    ind = closest_line(info, xy, ri, zi)

    info = info[ind]
  ENDIF ELSE info = info[0]
  
  CONTOUR, grid.psi, grid.r, grid.z, nlevels=50, xtit="R [m]" , ytit="Z [m]", /iso, color=1
  oplot_contour, info, xy, grid.r[*,0], grid.z[0,*], /periodic, color=2, thick=1.5
  
  ; Need to create boundary points uniformly spaced in theta
  rbndry = FLTARR(mfm)
  zbndry = FLTARR(mfm)
  
  FOR i=0, mfm-1 DO BEGIN
    theta = !DPI*2. * DOUBLE(i) / DOUBLE(mfm)
    
    ; Create a line going from magnetic axis outwards at angle theta
    r1 = [grid.rmagx, 100.*COS(theta)]
    z1 = [grid.zmagx, 100.*SIN(theta)]
    
    ; Get the location of the boundary in meters
    ri = REFORM(xy[0,info.offset:(info.offset+info.n-1)])
    zi = REFORM(xy[1,info.offset:(info.offset+info.n-1)])
    r2 = INTERPOLATE(grid.r[*,0], ri)
    z2 = INTERPOLATE(grid.z[0,*], zi)
    
    ; find where the radial line intersects the boundary
    hit = line_crossings(r1,z1, 0, r2, z2, 1, ncross=ncross)
    IF ncross NE 1 THEN BEGIN
      PRINT, "ERROR finding boundary"
      STOP
    ENDIF
    rbndry[i] = hit[0,0]
    zbndry[i] = hit[1,0]
    
    oplot, [grid.rmagx, rbndry[i]], [grid.zmagx, zbndry[i]], color=4
  ENDFOR
  
  oplot, rbndry, zbndry, psym=1
  
  OPENW, fid, "plasma_boundary.in", /GET_LUN
  FOR i=0, mfm-1 DO BEGIN
    PRINTF, fid, "  "+STR(rbndry[i], format="(E24.16)")+"     "+STR(zbndry[i], format="(E24.16)")
  ENDFOR
  CLOSE, fid

  minor_radius = MEAN(SQRT((rbndry - grid.rmagx)^2 + (zbndry - grid.zmagx)^2))
  eps = minor_radius / grid.rmagx ; inverse aspect ratio

  npts = N_ELEMENTS(grid.fpol)
  pn = FINDGEN(npts) / FLOAT(npts-1)

  ; Calculate Fprime and Pprime, taking derivative w.r.t. normalised psi
  pprime = DERIV(pn, grid.pres)
  fprime = DERIV(pn, grid.fpol)

  ffprime = grid.fpol * fprime

  ; Save pressure gradient
  OPENW, fid, "dp.in"
  FOR i=0, npts-1 DO BEGIN
    PRINTF, fid, STR(pprime[i])
  ENDFOR
  CLOSE, fid

  ; Save FFprime
  OPENW, fid, "fdf.in"
  FOR i=0, npts-1 DO BEGIN
    PRINTF, fid, STR(ffprime[i])
  ENDFOR
  CLOSE, fid
  
  OPENW, fid, "q_in.in"
  FOR i=0, npts-1 DO BEGIN
    PRINTF, fid, STR(grid.qpsi[i])
  ENDFOR
  CLOSE, fid

  ; Save XML input file
  
  IF KEYWORD_SET(new) THEN BEGIN
    
    OPENW, fid, "input_helena.xml"
    PRINTF, fid, "<?xml version='1.0'?>"
    PRINTF, fid, "<?xml-stylesheet type='text/xsl' href='./input_helena.xsl'"
    PRINTF, fid, "charset='ISO-8859-1'?>"
    PRINTF, fid, "<parameters>"
    
    PRINTF, fid, "<!-- shape parameters -->"
    PRINTF, fid, "<parameter>"
    PRINTF, fid, "  <name>ishape</name>"
    PRINTF, fid, "  <value>2</value>"
    PRINTF, fid, "  <description>Read plasma boundary from plasma_boundary.in</description>"
    PRINTF, fid, "</parameter>"
    PRINTF, fid, "<parameter>"
    PRINTF, fid, "  <name>mfm</name>"
    PRINTF, fid, "  <value>"+STR(mfm)+"</value>"
    PRINTF, fid, "  <description>Number of boundary points in theta</description>"
    PRINTF, fid, "</parameter>"
    PRINTF, fid, "<parameter>"
    PRINTF, fid, "  <name>isol</name>"
    PRINTF, fid, "  <value>0</value>"
    PRINTF, fid, "  <description>do NOT use Soloviev equilibrium</description>"
    PRINTF, fid, "</parameter>"
    PRINTF, fid, "<parameter>"
    PRINTF, fid, "  <name>ias</name>"
    PRINTF, fid, "  <value>1</value>"
    PRINTF, fid, "  <description>equilibrium is asymmetric</description>"
    PRINTF, fid, "</parameter>"
    PRINTF, fid, "<parameter>"
    PRINTF, fid, "  <name>npts</name>"
    PRINTF, fid, "  <value>"+STR(npts)+"</value>"
    PRINTF, fid, "  <description>Number of pprime and ffprime points in psi</description>"
    PRINTF, fid, "</parameter>"

    PRINTF, fid, "<!-- profile parameters -->"
    PRINTF, fid, "<parameter>"
    PRINTF, fid, "  <name>p%type</name>"
    PRINTF, fid, "  <value>7</value>"
    PRINTF, fid, "  <description>Read p-prime(psi) from dp.in</description>"
    PRINTF, fid, "</parameter>"
    PRINTF, fid, "<parameter>"
    PRINTF, fid, "  <name>gam%type</name>"
    PRINTF, fid, "  <value>7</value>"
    PRINTF, fid, "  <description>Read f * f-prime(psi) from fdf.in</description>"
    PRINTF, fid, "</parameter>"
    PRINTF, fid, "<parameter>"
    PRINTF, fid, "  <name>cur%type</name>"
    PRINTF, fid, "  <value>0</value>"
    PRINTF, fid, "  <description>fdf instead of current profile</description>"
    PRINTF, fid, "</parameter>"
    
    PRINTF, fid, "<!-- global parameters -->"
    PRINTF, fid, "<parameter>"
    PRINTF, fid, "  <name>eps</name>"
    PRINTF, fid, "  <value>"+STR(eps)+"</value>"
    PRINTF, fid, "  <description>Inverse aspect ratio</description>"
    PRINTF, fid, "</parameter>"
    PRINTF, fid, "<parameter>"
    PRINTF, fid, "  <name>alfa</name>"
    PRINTF, fid, "  <value>"+STR(psin)+"</value>"
    PRINTF, fid, "  <description>Normalised psi</description>"
    PRINTF, fid, "</parameter>"
    

    PRINTF, fid, "<!-- numerical parameters -->"
    PRINTF, fid, "<parameter>"
    PRINTF, fid, "  <name>nr</name>"
    PRINTF, fid, "  <value>101</value>"
    PRINTF, fid, "  <description>number of radial grid points</description>"
    PRINTF, fid, "</parameter>"
    PRINTF, fid, "<parameter>"
    PRINTF, fid, "  <name>np</name>"
    PRINTF, fid, "  <value>129</value>"
    PRINTF, fid, "  <description>number of poloidal grid points = 2^(n_1+1) for asymmetric equilibria</description>"
    PRINTF, fid, "</parameter>"
    PRINTF, fid, "<parameter>"
    PRINTF, fid, "  <name>nrmap</name>"
    PRINTF, fid, "  <value>"+STR(nrmap)+"</value>"
    PRINTF, fid, "  <description>number of radial points for mapping to final grid</description>"
    PRINTF, fid, "</parameter>"
    PRINTF, fid, "<parameter>"
    PRINTF, fid, "  <name>npmap</name>"
    PRINTF, fid, "  <value>"+STR(npmap)+"</value>"
    PRINTF, fid, "  <description>number of poloidal points for mapping to final grid. Must be 2^(n_2+1) for asymmetric equilibria</description>"
    PRINTF, fid, "</parameter>"
    PRINTF, fid, "<parameter>"
    PRINTF, fid, "  <name>nchi</name>"
    PRINTF, fid, "  <value>"+STR(npmap-1)+"</value>"
    PRINTF, fid, "  <description>Number of poloidal points written to equilibrium file. 2^n_3 for asymmetric equilibria</description>"
    PRINTF, fid, "</parameter>"
    PRINTF, fid, "<parameter>"
    PRINTF, fid, "  <name>niter</name>"
    PRINTF, fid, "  <value>20</value>"
    PRINTF, fid, "  <description>Maximum number of iterations with fixed grid</description>"
    PRINTF, fid, "</parameter>"
    PRINTF, fid, "<parameter>"
    PRINTF, fid, "  <name>nmesh</name>"
    PRINTF, fid, "  <value>3</value>"
    PRINTF, fid, "  <description>Maximum number of iterations for current or q profile</description>"
    PRINTF, fid, "</parameter>"
    PRINTF, fid, "<parameter>"
    PRINTF, fid, "  <name>errcur</name>"
    PRINTF, fid, "  <value>1.0E-07</value>"
    PRINTF, fid, "  <description>Convergence criterion for ff'</description>"
    PRINTF, fid, "</parameter>"
    
    PRINTF, fid, "<!-- diagnostic parameters -->"

    PRINTF, fid, "<parameter>"
    PRINTF, fid, "  <name>npr1</name>"
    PRINTF, fid, "  <value>1</value>"
    PRINTF, fid, "  <description>Output profiles and grid to out_he and SI output</description>"
    PRINTF, fid, "</parameter>"
    PRINTF, fid, "<parameter>"
    PRINTF, fid, "  <name>npr2</name>"
    PRINTF, fid, "  <value>1</value>"
    PRINTF, fid, "  <description>Output EQDSK file</description>"
    PRINTF, fid, "</parameter>"
    PRINTF, fid, "<parameter>"
    PRINTF, fid, "  <name>nprxmgr</name>"
    PRINTF, fid, "  <value>1</value>"
    PRINTF, fid, "  <description>Controls output of profiles to xmgr-readable ASCII files</description>"
    PRINTF, fid, "</parameter>"

    PRINTF, fid, "</parameters>"
    CLOSE, fid
  ENDIF ELSE BEGIN
    
    ; Old-style XML input
    OPENW, fid, "input_helena.xml"
    
    PRINTF, fid, '<?xml version="1.0"?>'
    PRINTF, fid, '<?xml-stylesheet type="text/xsl" href="./input_helena.xsl"'
    PRINTF, fid, 'charset="ISO-8859-1"?>'
    
    PRINTF, fid, "<parameters>"

    PRINTF, fid, "<!-- profile parameters -->"
    PRINTF, fid, " <profile_parameters>"
    PRINTF, fid, "   <p>"
    PRINTF, fid, "     <type> 7 </type>"
    PRINTF, fid, "   </p>"
    PRINTF, fid, "   <gam>"
    PRINTF, fid, "     <type> 7 </type>"
    PRINTF, fid, "   </gam>"
    PRINTF, fid, "   <cur>"
    PRINTF, fid, "     <type> 0 </type>"
    PRINTF, fid, "   </cur>"
    PRINTF, fid, "   <npts> "+STR(npts)+" </npts>"
    PRINTF, fid, " </profile_parameters>"
    
    PRINTF, fid, " <shape_parameters>"
    PRINTF, fid, "   <ishape> 2 </ishape>"
    PRINTF, fid, "   <isol> 0 </isol>"
    PRINTF, fid, "   <ias> 1 </ias>"
    PRINTF, fid, "   <mfm> "+STR(mfm)+" </mfm>"
    PRINTF, fid, " </shape_parameters>"
    
    PRINTF, fid, "<!-- global parameters -->"
    PRINTF, fid, " <global_parameters>"
    PRINTF, fid, "   <eps> "+STR(eps)+" </eps>"
    PRINTF, fid, "   <alfa> "+STR(psin)+" </alfa>"
    PRINTF, fid, " </global_parameters>"
    
    PRINTF, fid, "<!-- numerical parameters -->"
    PRINTF, fid, "  <numerical_parameters>"
    PRINTF, fid, "    <nr> 101 </nr>"
    PRINTF, fid, "    <np> 129 </np>"
    PRINTF, fid, "    <nrmap> "+STR(nrmap)+" </nrmap>"
    PRINTF, fid, "    <npmap> "+STR(npmap)+" </npmap>"
    PRINTF, fid, "    <nchi> "+STR(npmap-1)+" </nchi>"
    PRINTF, fid, "    <niter> 30 </niter>"
    PRINTF, fid, "    <nmesh> 5 </nmesh>"
    PRINTF, fid, "    <errcur> 0.1000E-06 </errcur>"
    PRINTF, fid, "  </numerical_parameters>"

    PRINTF, fid, "<!-- diagnostics parameters -->"
    PRINTF, fid, " <diagnostics_parameters>"
    PRINTF, fid, "   <npr1> 1 </npr1>"
    PRINTF, fid, "   <npr2> 1 </npr2>"
    PRINTF, fid, "   <nprxmgr> 1 </nprxmgr>"
    PRINTF, fid, "   <ndiag> 1 </ndiag>"
    PRINTF, fid, "   <npl1> 0 </npl1>"
    PRINTF, fid, " </diagnostics_parameters>"
    
    PRINTF, fid, "</parameters>"
    CLOSE, fid
  ENDELSE
  
  PRINT, "DONE"
  STOP
END

