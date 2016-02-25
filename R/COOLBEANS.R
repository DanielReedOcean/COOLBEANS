###############################################################################
# THIS CODE IS DISTRIBUTED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OR 
# CONDITIONS OF ANY KIND - WITHOUT EVEN THE IMPLIED WARRANTY OF MERCHANTABILITY
# OR FITNESS FOR A PARTICULAR PURPOSE.
###############################################################################


###############################################################################
# Purpose: Model of sediment oxygen demand & auxilliary functions
# Author:  Daniel Reed (dan.reed@wsu.edu)
# Date:    January, 13th 2016
###############################################################################


# Arguments: water depth (m)
# Returns: organic carbon burial rate from Middelburg et al. (1997) 
#(mol/[m2 y])
C.burial <- function(z) {
  return(10 * 4.4 * 10^{-0.84672973 - 0.00061506 * z})
}

# Arguments: water depth (m)
# Returns: organic carbon remineralisation rate from Middelburg et al. (1997) 
# (mol/[m2 y])
C.remin <- function(z) {
  return(10 * 1.8 * 10^{-0.50860503 - 0.000389 * z})
}

# Arguments: depth (m)
# Returns: proportion of organic flux that is remineralised (unitless)
p.remin <- function(z) {
  return(C.remin(z) / (C.remin(z) + C.burial(z)))
}

# Arguments: net primary production (mol/[m2 y]); water depth (m)
# Returns: export flux attenuated by depth relationship from 
# Andersson et al. (2004) to give organic carbon flux as a function of 
# depth for a particular NPP (mol/[m2 y])
C.flux <- function(NPP, z) {
  return( NPP * ( (1-0.17) * exp(-0.018*z) + 0.17 * exp(-0.00046*z) ) )
}

# Arguments: water depth (m); net primary production (mol/[m2 y]); 
# water column supported denitrification (mol/[m2 y])
# Returns: sediment oxygen demand model (mol/[m2 y])
BOD <- function(z, NPP, C.flux.extra = 0) {
  # Calculate organic carbon flux to sediments from net primary production (NPP)
  # and water depth
  C.flux.sed <- C.flux(NPP, z) + C.flux.extra
    
  # Calculate benthic oxygen demand
  Oxygen.demand <- (((77/60) * p.remin(z) - (17/60)) * C.flux.sed)  
  
  #No export, no oxygen uptake
  if(any(NPP <= 0)) warning(paste("Net primary production is ", NPP[NPP <= 0],
                                  " - must be positive.\n", sep=""))
  
  if(any(z <= 0)) warning(paste("Water depth is ", z[z <= 0],
                                " - must be positive.\n", sep=""))
  
  # Return benthic oxygen demand
  return(Oxygen.demand)
}

# Arguments: benthic oxygen demand (mol/[m2 y]); surface oxygen (mol/m3);  
# bottom water oxygen (mol/m3); bottom layer thickness (m)
# Returns: exchange coefficient (/y)
Exchange.coefficient <- function(BOD, O2.surf, O2.bottom, L = 1){
  # Calculate exchange coefficient
  Exchange <- BOD/((O2.surf - O2.bottom) * L)
  
  # Mask array for those values that are undefined
  mask <- (0.05 > (O2.surf - O2.bottom)/O2.surf | O2.bottom > O2.surf)
  
  # Update undefined values using expression from scaling analysis
  Exchange[mask] <- BOD[mask]/(0.1 * O2.surf[mask] * L)
  
  # Return exchange coefficients
  return(Exchange)
}

# Arguments: benthic oxygen demand (mol/[m2 y]); surface oxygen (mol/m3);
# exchange coefficient (/y); bottom layer thickness (m)
# Returns: bottom water oxygen (mol/m3)
Bottom.water.O2 <- function(O2.surf, BOD, Exchange, L = 1){
  # Calculate bottom water oxygen
  O2.bottom <- O2.surf - BOD/(Exchange * L)
  
  # Return bottom water oxygen
  return(O2.bottom)
}

# Arguments: benthic oxygen demand (mol/[m2 y]); exchange coefficient (/y) 
# bottom layer thickness (m)
# Returns: Sensitivity to exchange coefficient (-)
Exchange.sensitivity <- function(BOD, Exchange, L = 1){
  # Return sensitivity of bottom water oxygen to vertical mixing
  return(BOD/(L * Exchange^2))
}

# Arguments: water depth (m); exchange coefficient (/y)  
# Returns: Sensitivity to terrestrial nutrient loading (-)
Loading.sensitivity <- function(Depth, Exchange, L = 1){
  # Calculate proportion of organic matter remineralisaed
  gamma <- p.remin(Depth)
  
  # Return sensitivity of bottom water oxygen to nutrient loading
  return(-10 * ((77/60) * gamma - (17/60))/
           (Exchange * L * ((733/1500) * gamma + (1007/1500))))
}

# Arguments: net primary production (mol C/[m2 y]); exchange coefficient (/y);
# surface water oxygen (mol/m3); surface water temperature (ºC); water depth (m);
# area of COSCAT (m2); bottom layer thickness (m)
# Returns: Maximum hypoxic area (m2)
A.hy <- function(NPP, Exchange, O2.surf, Depth, Area, L = 1){
  # Calculate total depositional flux of organic matter in COSCAT
  M.sed <- C.flux(NPP, Depth) * Area
  
  # Calculate proportion of organic matter remineralisaed
  gamma <- p.remin(Depth)
  
  # Return maximum area of hypoxia supported by depositional flux
  return(M.sed * (gamma * 77/60 - 17/60) / (L * Exchange * (O2.surf - 0.063)))
}
