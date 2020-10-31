;;; Aggregates 3-hourly NARCCAP data to averages over daily and longer periods

;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"



;; This stuff is localized into functions so that seasons or their
;; ordering can be redefined.  Functions could be rewritten to take an
;; optional argument for whether you wanted monsoons, etc.  May move
;; into a separate file at some point.

; Example: print(days_in_year(1996, "julian")) ; 366
function days_in_year(year, calendar)
begin
  z0 = new(12,"integer")
  z0 = 0
  years = z0 + year
  years@calendar = calendar
  months = ispan(1,12,1)  
  return( sum(days_in_month(years,months)) ) 
end


function getnames(interval, calendar)
begin
  if (interval .eq. "season") then
    return((/"error","DJF","MAM","JJA","SON"/))
  end if
  if (interval .eq. "month") then
    return((/"error","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/))
  end if
  if (interval .eq. "day") then
      return( array_append_record("error", sprintf("%03.0f", ispan(1,days_in_year(1996,calendar),1)),0 ) )
  end if 
  return("getnames() error: unrecognized interval")
end


function mon2seas(month:numeric)
local season
begin
  season = new(dimsizes(month),"integer")
  season = where(month.eq.12 .or. month.eq. 1 .or. month.eq. 2, 1, season)  ;; DJF
  season = where(month.eq. 3 .or. month.eq. 4 .or. month.eq. 5, 2, season)  ;; MAM
  season = where(month.eq. 6 .or. month.eq. 7 .or. month.eq. 8, 3, season)  ;; JJA
  season = where(month.eq. 9 .or. month.eq.10 .or. month.eq.11, 4, season)  ;; SON
  return(season)
end


;;  For converting times to dates

format = "%Y-%N-%D %H:%M:%S"


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;; Command-line inputs:
;; infile
;; outfile
;; varname
;; interval: day, month, season, or year (no default)
;; method:  mean, min, max (defaults to mean)
;; outtime: start, mid, end  (defaults to end)
;; cyclic  (boolean option: average over multiple years?)
;; check   (boolean option: print provenance checking info?)
;; offset  (optional: float added to time to shift day boundary)
;; taint   (boolean option: propagate missing values?)
;; verbose (boolean option: print progress indicators?)
;; donemsg (boolean option: print completion message?)
;; drop29  (boolean option: drop Feb 29th in leapyears, this is the default)
;; drop366 (boolean option: drop day 366 in leapyears)

if (.not. isvar("drop366")) then
  drop29 = True
end if

if (.not. isvar("verbose")) then
  verbose = False
end if

if (.not. isvar("donemsg")) then
  donemsg = False
end if

if (.not. isvar("check")) then
  check = False
end if

if (.not. isvar("cyclic")) then
  cyclic = False
end if

if (.not. isvar("taint")) then
  taint = False
end if

if (.not.isvar("outtime")) then
  outtime = "end"
end if
if (outtime.ne."start" .and. outtime.ne."mid" .and. outtime.ne."end") then
  print("Invalid outtime: '"+outtime+"'.  Use start, mid, or end.")
  exit
end if


if (.not.isvar("method")) then
  method = "mean"
end if

if (method.ne."mean" .and. method.ne."min" .and. method.ne."max") then
  print("Invalid method: '"+method+"'.  Use mean, min, or max.")
  exit
end if

if (method .eq. "max") then
  mstring = "maximum"
else if (method .eq. "min") then
    mstring = "minimum"
  else
    mstring = method
  end if
end if

if (.not.isvar("interval")) then
  print("Specify interval: day, month, season, or year")
  exit
end if
if (interval.ne."day" .and. interval.ne."month" .and. interval.ne."season" .and. interval.ne."year") then
  print("Invalid interval: '"+interval+"'.  Use day, month, season, or year.")
  exit
end if

if(cyclic .and. interval .eq. "year") then
  print("Error: annual interval not allowed with cyclic option.")
  ;; actually, this is just averaging over all time.  The only reason to
  ;; disallow it is clarity.  If we add an "all" option, we can make
  ;; this case synonymous with that and implement it with just a
  ;; dim_avg_n shortcut.
  exit
end if

;if(cyclic .and. interval .eq. "day") then
;  print("Error: day interval not allowed with cyclic option (due to leap days).")
  ;; once NCL v.6.1.0 is readily available, this is actually an
  ;; allowable option.  Use day_of_year() to convert from Y/M/D output
  ;; of cd_calendar to julian day, and just discard day 366 in leap
  ;; years.  (day_of_year doesn't recognize calendars before v.6.1.0)
;  exit
;end if

if (verbose) then print("Opening "+infile) end if

fin = addfile(infile, "r")

if (verbose) then print("Reading time coordinate") end if

time = fin->time
nt = dimsizes(time)

if(isvar("offset")) then 
  time = time + offset
end if

if(.not. any(typeof(time).eq.(/"float","double"/))) then
  print("Error: time is not floating-point.")
  exit
end if


;; get or create time@bounds as appropriate

if (verbose) then print("Setting up time_bnds") end if

if (isatt(time,"bounds")) then
  tbname = time@bounds
  time_bnds = fin->$tbname$
else
  time_bnds = transpose(onedtond(time,(/2,nt/)))
end if

 
;; Extract month & year using ut_calendar

utc = cd_calendar(time, 0)

year  = floattointeger(utc(:,0))
month = floattointeger(utc(:,1))
day   = floattointeger(utc(:,2))
year@calendar = fin->time@calendar
doy   = day_of_year(year,month,day)

; remove feb 29 from doy
; note that isleapyear() automatically resolves different calendars, given the calendar attribute
; https://www.ncl.ucar.edu/Document/Functions/Built-in/isleapyear.shtml
isleap = isleapyear(year)
if (drop29) then
  do i = 0, nt-1
    if( isleap(i) .and. doy(i).ge.59 ) then
      doy(i) = doy(i)-1
    end if
  end do
end if

if(interval .eq. "day") then
  tunits = str_get_field(time@units,1," ")
  if (tunits .eq. "days") then
;;      jday  = doubletointeger(time)
      jday = floor(time)
    else if (tunits .eq. "hours") then
;;      jday  = doubletointeger(time/24)
      jday = floor(time/24)
    else if (tunits .eq. "minutes") then
        jday = floor(time/60/24)
    else
      print("Error: don't know how to convert time unit '"+tunits+"'into julian day.")
      exit
    end if
    end if
  end if
end if


season = mon2seas(month)
names =  getnames(interval, fin->time@calendar)


;; Construct mapping info:

if (verbose) then print("Constructing index mapping") end if

;; oci == output cell index array

oci = new(nt,"integer")

if (cyclic) then
  ;; climatology - multi-year averages

  if(interval .eq. "day") then
    nout = dimsizes(names) - 1
    bin = doy
  end if

  if(interval .eq. "month") then 
    nout = dimsizes(names) - 1
    bin = month
  end if

  if (interval .eq. "season") then
    nout = dimsizes(names) - 1
    bin = season   
  end if

  ;; The -1 is because January = month 1 goes in array position time(0), etc.
  do i = 1, nout
    oci(ind(bin .eq. i)) = i-1
  end do

else
  ;; normal (linear/non-cyclic) intervals
  
 if(interval .eq. "day"   ) then bin = jday   end if
 if(interval .eq. "month" ) then bin = month  end if
 if(interval .eq. "season") then bin = season end if
 if(interval .eq. "year"  ) then bin = year   end if
 
 nout = 0
 oci(0) = 0 
 do i = 1, nt-1
   if (bin(i) .ne. bin(i-1)) then
     nout = nout + 1
   end if
   oci(i) = nout
 end do
 nout = nout + 1
end if



;; And now for the actual averaging of data

if (verbose) then print("Setting up result array") end if

;; This is a cheap way to get metadata from a non-copied var
outdata = fin->$varname$(0:nout-1,:,:)
outdata(:,:,:) = outdata@_FillValue
odim = dimsizes(outdata)

misscount = new(nout,"integer")
stepcount = new(nout,"integer")
cell_bnds = new((/nout,2/), "double","No_FillValue") ;; will become time_bnds; follows CF rules
date = new (nout, "string")                          ;; plain-text version of output time coordinate


;; so the deletes don't error...
map  = 0
slab = 0 

;; Note that we dodge around reading in fin->$varname$ so that we
;; don't crash in cases where it's too big to fit into memory.

if (verbose) then print("Aggregating") end if

do i=0, nout-1

  delete(map)
  map = ind(oci .eq. i)

  if(dimsizes(map) .eq. 1) then

    ;; A little kludgey.  Probably ought to do something with expected
    ;; of timesteps and a threshold for allowable %missing instead of
    ;; simple yes-no tainting by missing value...

    if(.not. taint) then
      outdata(i,:,:) = fin->$varname$(map(0):map(0),:,:)
    end if
    misscount(i) = num(ismissing(outdata(i,:,:)))
    stepcount(i) = 1
    cell_bnds(i,0) = time_bnds(map(0),0)
    cell_bnds(i,1) = time_bnds(map(0),1)

  else 

    delete(slab)
    slab = fin->$varname$(map,:,:)
  
    if(method.eq."mean") then
      outdata(i,:,:) = dim_avg_n(slab,0)
    end if
    if(method.eq."min") then
      outdata(i,:,:) = dim_min_n(slab,0)
    end if
    if(method.eq."max") then
      outdata(i,:,:) = dim_max_n(slab,0)
    end if
  
    if(taint) then
      ds = dimsizes(slab)
      do tt = 0,ds(0)-1
        outdata(i,:,:) = outdata(i,:,:) + (slab(tt,:,:)*0)
      end do
      delete(ds)
    end if 
  
    nmap = dimsizes(map)
    misscount(i) = num(ismissing(slab))
    stepcount(i) = nmap
    cell_bnds(i,0) = time_bnds(map(0),0)
    cell_bnds(i,1) = time_bnds(map(nmap-1),1)

  end if

  if(cyclic) then
    date(i) = names(i+1)
  else
    if (interval .eq. "day") then
      date(i) = ""+year(map(0))+"-"+month(map(0))+"-"+day(map(0))
    else 
      if (interval .eq. "year") then
        date(i) = year(map(0))
      else
        date(i) = names(bin(map(0)))+" "+year(map(0))
      end if
    end if
  end if
end do

if (verbose) then print("Aggregating") end if


if(isvar("offset")) then 
  cell_bnds = cell_bnds - offset
end if

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; Done with the actual aggregation.  From here on out, we're just
;; generating output

;; Add appropriate metadata to output variables

if (verbose) then print("Adding metadata to results") end if

delete(time)
delete(time_bnds)

if(outtime.eq."start") then
    time = cell_bnds(:,0)
end if
if(outtime.eq."mid") then
  if(cyclic) then
    outyear = tointeger((year(0)+year(nt-1))/2)
    if(interval.eq."month") then
      z0 = new(12,"integer")
      z0 = 0
      z0@calendar = fin->time@calendar
      mout = ispan(1,12,1)
      time = cd_inv_calendar(z0+outyear,mout,z0+15,z0+12,z0,z0,fin->time@units, z0(0))
    end if
    if(interval.eq."season") then
      z0 = new(4,"integer")
      z0 = 0
      z0@calendar = fin->time@calendar
      mout = (/1,4,7,10/)
      time = cd_inv_calendar(z0+outyear,mout,z0+15,z0+12,z0,z0,fin->time@units, z0(0))
    end if
    if(interval.eq."day") then
      maxdays = days_in_year(1996,fin->time@calendar)
      z0 = new(maxdays,"integer")
      z0 = 0
      z0@calendar = fin->time@calendar
      time = cd_inv_calendar(z0+outyear,z0+1,ispan(1,maxdays,1),z0+12,z0,z0,fin->time@units, z0(0))
    end if
  else
    time = cell_bnds(:,1)
    time = (time + cell_bnds(:,0))/2
  end if
end if
if(outtime.eq."end") then
  time = cell_bnds(:,1)
end if
if (cyclic) then
time@climatology = "time_bnds"
else
time@bounds = "time_bnds"
end if
time@axis = "T"
time@calendar = fin->time@calendar
time@units = fin->time@units
time@standard_name = "time"
time@long_name = "time"
;;delete(time@_FillValue)
time!0 = "time"
time&time = time

time_bnds = cell_bnds
time_bnds!0 = "time"
time_bnds!1 = "bnds"
;;delete(time_bnds@_FillValue)


if (cyclic) then
  cmstring = "time: "+mstring+" within years time: "+mstring+" over years"
else
  cmstring = "time: "+mstring
end if
if (isatt(outdata, "cell_methods")) then
  cmstring = outdata@cell_methods + " " + cmstring
end if
outdata@cell_methods = cmstring


;; Create output file

if (verbose) then print("Creating output file") end if

system("rm -f "+outfile)
fout = addfile(outfile, "c")
filedimdef(fout,"time",-1,True)  ;; make time dimension unlimited


if (verbose) then print("Copying other variables from input") end if

;; copy global attributes

att_names = getvaratts(fin) 

do i = 0,dimsizes(att_names)-1
  fout@$att_names(i)$ = fin@$att_names(i)$
end do

;; copy over other variables

var_names = getfilevarnames (fin) ;

do i = 0,dimsizes(var_names)-1      

  if (var_names(i) .eq. "time")      then continue end if
  if (var_names(i) .eq. "time_bnds") then continue end if
  if (var_names(i) .eq. varname)     then continue end if
  if (var_names(i) .eq. "date")      then continue end if
  if (var_names(i) .eq. "date_bnds") then continue end if

  fout->$var_names(i)$ = fin->$var_names(i)$
end do


if (verbose) then print("Writing results to output") end if

;; add our data

fout->time_bnds = time_bnds
fout->$varname$ = outdata
fout->time      = time
; For some reason, time needs to come *last* or it gets values from fin instead
delete(fout->time@_FillValue)  ;; coming from filedimdf call to make it unlimited, I think


;; append history entry

hstring = systemfunc("date")
hstring = hstring + ": ncl aggregate.ncl infile="+infile+" outfile="+outfile+" varname="+varname+" interval="+interval
if (cyclic) then
  hstring = hstring + " cyclic=True"
end if 
;;if (midpoint) then
;;  hstring = hstring + " midpoint=True"
;;end if 
hstring = hstring+" method="+method
hstring = hstring+" outtime="+outtime
fout@history = hstring+inttochar(10)+fout@history

;; update unique tracking_id
fout@tracking_id = systemfunc("uuidgen")


if (donemsg) then print("Done: "+outfile) end if

;; provenance checking - to be sure that your data is coming from where you think it should.

if(check) then

if (verbose) then print("") end if
if (verbose) then print("") end if
if (verbose) then print("Provenance check:") end if
if (verbose) then print("") end if

  db = ndtooned(time_bnds)
  db@calendar = time@calendar
  db@units = time@units
  date_bnds = onedtond(cd_string(db,format), dimsizes(time_bnds))
  
  tab = "	"

  if(cyclic) then spacer = "" else spacer = "        " end if

  print("#index"+tab+"#date"+spacer+tab+"#starting_cell_bdy"+tab+"#ending_cell_bdy"+tab+"#nsteps"+tab+"#missing")

  nxy = product(dimsizes(outdata(0,:,:)))

  do i=0,nout-1
    print(i+tab+date(i)+tab+date_bnds(i,0)+tab+date_bnds(i,1)+tab+stepcount(i)+tab+misscount(i)/nxy)
  end do

  print("")
  print("")

  if(cyclic) then

    nb = dimsizes(names)

    head = ""
    do b = 1, nb-1
      head = head+tab+names(b)
    end do
    print(""+head)

    do y = min(year), max(year)
      line = ""+y
      do b = 1, nb-1
        line = line + tab + num(year.eq.y .and. bin.eq.b)
      end do
      print(""+line)
    end do

  end if


end if

exit


;; Copyright 2009-2012 Univ. Corp. for Atmos. Research
;; Author: Seth McGinnis, mcginnis@ucar.edu
