# Function to calculate saturated vapor pressure in kPa
# Requires Tmax, Tmin in degrees C
# Tmax - max temp; Tmin - min temp, 
svp <- function(Tmax,Tmin)
{
  esmax <- 0.6108*exp(17.27 * Tmax/(Tmax + 237.3))
  esmin <- 0.6108*exp(17.27 * Tmin/(Tmin + 237.3))
  et <- (esmax + esmin)/2
  return(et)
}

# Function to calculate actual vapor pressure with dew point
# Requires dew point temperature in degree C; eadew in kPa
eadew <- function(Tdew)
{
  ead <- 0.6108*exp(17.27 * Tdew/(Tdew + 237.3))
  return(ead)
}

# Function to calculate vapor pressure deficit vpd
# es and ea in kPa
#  ea - actual vapor pressure and es - saturated vapor pressure

vpd <- function(es,ea)
{
  vpdx <- es-ea
  return(vpdx)
}

# Function to calculate actual vapor pressure with min temperature
# Requires Tmin - Min. Temperaure in degree C and AFLAG 
# AFLAG = 1 if arid or semi-arid and 0 if humid
eatmin <- function(Tmin,AFLAG)
{
  if(AFLAG == 1)
  {
    Tmin <- Tmin - 2.0
  }
  eat <- 0.6108*exp(17.27 * Tmin/(Tmin + 237.3))
  return(eat)
}

# Function to calculate actual vapor pressure with Av. Relative Humidity
# Requires Tmax, Tmin in degree C and RH in %
eaRHav <- function(RH,Tmax,Tmin)
{
  es <- svp(Tmax,Tmin)
  ear <- (RH/100)*svp
  return(ear)
}

# Function to calculate actual vapor pressure with min and max Relative Humidity
# Requires Tmin, Tmax, RHmin (%) and RHmax (%)
eaRHmm <- function(RHmin, RHmax, Tmax, Tmin)
{
  etmin <- 0.6108*exp(17.27 * Tmin/(Tmin + 237.3))
  etmax <- 0.6108*exp(17.27 * Tmax/(Tmax + 237.3))
  earmm <- 0.5 *((RHmin/100)*etmax + (RHmax/100) * etmin)
  return(earmm)
}


# Function to slope of the VP-Temperature curve (kPa/C)
# Requires Tmax and Tmin
deltaVP <- function(Tmax,Tmin)
{
  T <- (Tmax + Tmin)/2
  num <- 0.6108*exp(17.27*T/(T+237.3))
  dnum <- (T + 237.13)^2
  delta <- 4098*(num/dnum)
  return(delta)
}

# Function for Psychometric Constant (kPa/C)
# z - elevation in m above MSL
psychC <- function(z)
{
  P <- 101.3 * ((293-0.0065*z)/293)^5.26
  gamma <- 0.665 * 10^-3 * P
  return(gamma)
}

# Windspeed at 2 m height
# Requires u velocity (m/s) at height y (m) above ground surface
u2 <- function(u,y)
{
  u2 <- u * (4.87)/(log(67.8*y - 5.42))
  return(u2)
}

# Extraterrstrial Radiation MJ/m^2/d
# Requires, Latitude (decimal degrees); and Date (date object)
Ra <- function(lat,datex)
{
  J <- as.numeric(strftime(datex,format="%j"))
  dr <- 1 + 0.033*cos(2*pi*J/365)
  dc <- 0.409*sin(2*pi*J/365 - 1.39)
  # convert latitude to radians
  latr <- lat*pi/180
  
  # Sunset angle
  ws <- acos(-tan(latr)*tan(dc))
  
  # Calculate extraterrestrial solar radiation
  num1 <- ws*sin(latr)*sin(dc)+cos(latr)*cos(dc)*sin(ws)
  Gsc <- 0.0820
  Rax <- (24*(60)/pi)*Gsc*dr*num1
  return(Rax)
}

# Net Shortwave incoming solar Radiation
# Requires total solar radiation (MJ/m^2/d) and albedo (default = 0.23)
Rns <- function(Rs,alpha=0.23)
{
  Rnsx <- (1-alpha)*Rs
  return(Rnsx)
}



# Net Longwave solar Radiation
# Compute net outgoing long-wave radiation (Rnl)
Rnl <- function(Tmax,Tmin,ea,z,Ra,Rs)
{
  # Calculate clear-sky radiation Rso
  Rso = (0.75 + 2E-5*z)*Ra
  sigma <- 4.903 * 10^-09
  num1 <- ((Tmax+273.15)^4 + (Tmin+273.15)^4)/2
  num2 <- (0.34-0.14*sqrt(ea))
  num3 <- 1.35*(Rs/Rso)-0.35
  Rnlx <- sigma*num1*num2*num3
  return(Rnlx)
}

# Compute Net Radiation Rn
# Net Incoming SHort-wave radiation minus Outgoing longwave radiation
Rn <- function(Tmax,Tmin,ea,z,Ra,Rs,alpha=0.23)
{
  Rnsx <- Rns(Rs,alpha)
  Rnlx <- Rnl(Tmax,Tmin,ea,z,Ra,Rs)
  Rnx <- Rnsx - Rnlx
  return(Rnx)
}

# Function to Compute Penman-Monteith Equation
# Requires separate calculation of extra-terrestrial solar radiation
# Requires separate calculation of actual vapor pressure

pm <- function(Tmax,Tmin,u2,Rs,Ra,ea,lat,z,alpha)
{
  T <- (Tmax+Tmin)/2
  DELTA <- deltaVP(Tmax,Tmin)
  gamma <- psychC(z)
  Rnx <-  Rn(Tmax,Tmin,ea,z,Ra,Rs,alpha=0.23)
  et <- svp(Tmax,Tmin)
  vpdx <- vpd(et,ea)
  DNUM <- DELTA + gamma*(1 + 0.34*u2)
  NUMR <- 0.408*DELTA*Rnx
  NUMW <- gamma*(900/(T+273))*u2*vpdx
  ETrad <- NUMR/DNUM
  ETaer <- NUMW/DNUM
  ETtot <- ETrad + ETaer
  etall <- cbind(ETrad,ETaer,ETtot)
  colnames(etall) <- c('ETrad','ETaer','ETtot')
  return(etall)
}

# Function to Calculate Solar Radiation based on Hargreaves-Samani Eqn
# Rs is in MJ/m^2/d
# Requires Tmax, Tmin in degrees C, Extraterrestrial Solar Radiation R and
# FLG whose value is 1 for inland and 0 for coastal locations 
Rshar <- function(Tmax,Tmin,Ra,FLG)
{
  if(FLG == 1)
  {
    Ks <- 0.16}else{
      Ks <- 0.19
    }
  Rsh <- Ks * sqrt((Tmax - Tmin))*Ra
  return(Rsh)
}





