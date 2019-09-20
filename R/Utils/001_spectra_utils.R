#====================================================================================== -

#HEADER ----

# Author: Jonas Anderegg, ETH Zürich
# Copyright (C) 2019  ETH Zürich, Jonas Anderegg (jonas.anderegg@usys.ethz.ch)

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

#====================================================================================== -


# MAIN FUNCTIONS -----

#read spectra from .csv files
f_spc_read_fromPy <- function(dir, Exp = "FP") {
  
  #read .csv files
  setwd(dir)
  fnames <- as.list(list.files(dir, pattern = ".csv", full.names = TRUE))
  ASD <- lapply(fnames, read.csv)
  #extract filenames which will be used in later steps to create spc_ID
  #replace "Cnp" by "Can" to match 2016 and 2017 nomenclature
  names(ASD) <- fnames %>% unlist() %>% basename() %>% gsub(x = ., "Cnp", "Can")
  #extract wavelengths and reflectance
  rflt <- sapply(ASD, "[[", "rflt")
  wvlt <- sapply(ASD, "[[", "wvlt")
  names <- sapply(fnames, basename)
  
  #transpose to wide format
  df <- cbind(wvlt[,1], rflt)
  df <- t(df)[-1,] %>% as.data.frame()
  
  #rename columns
  colnames(df) <- unique(paste("rflt", wvlt, sep = "_"))
  
  #extract info from spectra name
  pn0 <- strsplit(names, "_", fixed = TRUE) %>% unlist()
  Plot_ID <- grep(Exp, pn0, value = TRUE)
  replicate <- stringr::str_sub(names, -5, -5)
  meas_date <- stringr::str_sub(strptime(grep("201", names, value = TRUE), "%Y%m%d"), 1, 10) %>% as.Date()
  
  #create data frame
  df <- data.frame(Plot_ID, meas_date, replicate, df)
  
  #drop reference measurements
  #drop unused factor levels
  #convert to tibble for better handling
  df <- df[!grepl("Ref", df$Plot_ID),]
  df$Plot_ID <- df$Plot_ID %>% droplevels()
  df <- tibble::as_tibble(df)
  
  return(df)
  
}

#smooth spectra and calculate derivatives
f_spc_smooth <- function(data, p = 3, w = 11, m = 0) {
  
  # data: a data frame
  # p: polynomial order (prospectr::savitzkyGolay)
  # m: differentiation order (prospectr::savitzkyGolay)
  # w: window size (must be odd) (prospectr::savitzkyGolay)
  
  #select reflectance data
  spc_raw <- data %>% dplyr::select(contains("rflt_"))
  #smooth reflectance
  spc_smth <- prospectr::savitzkyGolay(spc_raw, p = p, w = w, m = m) %>% 
    tibble::as.tibble()
  
  #combine measurement information with smoothed spectra
  df_smth <- data %>% dplyr::select(-contains("rflt_")) %>% 
    dplyr::bind_cols(., spc_smth)
  
  return(df_smth)
  
}

#average spectra
f_spc_avg <- function(data) {
  
  df_avg <- data %>% 
    dplyr::select(-contains("rflt_"), everything(), -replicate) %>%  
    group_by(Plot_ID, meas_date) %>% 
    dplyr::summarise_all(funs(mean)) %>% 
    ungroup()
  
} 

#compute spectral indices
f_calc_si <- function(data) {
  
  # BANDS ########################################################################## -
  
  ALI5 <- rowMeans(data[ ,match(paste("rflt_", 775:805, sep=""), colnames(data))])
  ALI6 <- rowMeans(data[ ,match(paste("rflt_", 845:890, sep=""), colnames(data))])
  ALI7 <- rowMeans(data[ ,match(paste("rflt_", 1200:1300, sep=""), colnames(data))])
  ASTER4 <- rowMeans(data[ ,match(paste("rflt_", 1600:1700, sep=""), colnames(data))])
  ASTER5 <- rowMeans(data[ ,match(paste("rflt_", 2145:2185, sep=""), colnames(data))])
  ASTER6 <- rowMeans(data[ ,match(paste("rflt_", 2185:2225, sep=""), colnames(data))])
  ASTER7 <- rowMeans(data[ ,match(paste("rflt_", 2235:2285, sep=""), colnames(data))])
  ASTER8 <- rowMeans(data[ ,match(paste("rflt_", 2295:2365, sep=""), colnames(data))])
  MODIS5 <- rowMeans(data[ ,match(paste("rflt_", 1230:1250, sep=""), colnames(data))])
  MODIS6 <- rowMeans(data[ ,match(paste("rflt_", 1628:1652, sep=""), colnames(data))])
  MODIS7 <- rowMeans(data[ ,match(paste("rflt_", 2105:2155, sep=""), colnames(data))])
  R1050_ASD <- rowMeans(data[ ,match(paste("rflt_", 1046:1054, sep=""), colnames(data))])
  R1090to1110 <- rowMeans(data[ ,match(paste("rflt_", 1090:1110, sep=""), colnames(data))])
  R1094_HYP <- rowMeans(data[ ,match(paste("rflt_", 1089:1099, sep=""), colnames(data))]) #B095 1094.09 10.99
  R1100_5 <- rowMeans(data[ ,match(paste("rflt_", 1098:1102, sep=""), colnames(data))]) #interpolated ASD
  R1180to1220 <- rowMeans(data[ ,match(paste("rflt_", 1180:1220, sep=""), colnames(data))])
  R1200_5 <- rowMeans(data[ ,match(paste("rflt_", 1198:1202, sep=""), colnames(data))]) #interpolated ASD
  R1205_HYP <- rowMeans(data[ ,match(paste("rflt_", 1201:1209, sep=""), colnames(data))]) #B106 1205.07 10.79
  R1220 <- rowMeans(data[ ,match(paste("rflt_", 1216:1224, sep=""), colnames(data))])
  R1225 <- rowMeans(data[ ,match(paste("rflt_", 1221:1229, sep=""), colnames(data))])
  R1240_AVIRIS <- rowMeans(data[ ,match(paste("rflt_", 1236:1244, sep=""), colnames(data))])
  R1250_ASD <- rowMeans(data[ ,match(paste("rflt_", 1246:1254, sep=""), colnames(data))])
  R1265to1285 <- rowMeans(data[ ,match(paste("rflt_", 1265:1285, sep=""), colnames(data))])
  R1300_100 <- rowMeans(data[ ,match(paste("rflt_", 1250:1350, sep=""), colnames(data))])
  R1450_100 <- rowMeans(data[ ,match(paste("rflt_", 1400:1500, sep=""), colnames(data))])
  R1510_AVIRIS <- rowMeans(data[ ,match(paste("rflt_", 1506:1514, sep=""), colnames(data))])
  R1555_ASD <- rowMeans(data[ ,match(paste("rflt_", 1551:1559, sep=""), colnames(data))])
  R1650 <- rowMeans(data[ ,match(paste("rflt_", 1649:1651, sep=""), colnames(data))])
  R1650_CS <- rowMeans(data[ ,match(paste("rflt_", 1550:1750, sep=""), colnames(data))])
  R1659_HYP <- rowMeans(data[ ,match(paste("rflt_", 1654:1664, sep=""), colnames(data))]) #B151 1659 11.53
  R1680_AVIRIS <- rowMeans(data[ ,match(paste("rflt_", 1676:1684, sep=""), colnames(data))])
  R1754_AVIRIS <- rowMeans(data[ ,match(paste("rflt_", 1750:1758, sep=""), colnames(data))])
  R2020_IRIS <- rowMeans(data[ ,match(paste("rflt_", 2019:2021, sep=""), colnames(data))])
  R2100_IRIS <- rowMeans(data[ ,match(paste("rflt_", 2099:2101, sep=""), colnames(data))])
  R2220_IRIS <- rowMeans(data[ ,match(paste("rflt_", 2219:2221, sep=""), colnames(data))])
  R360to370 <- rowMeans(data[ ,match(paste("rflt_", 360:370, sep=""), colnames(data))])
  R430_6 <- rowMeans(data[ ,match(paste("rflt_", 428:432, sep=""), colnames(data))])
  R440_SD <- rowMeans(data[ ,match(paste("rflt_", 438:442, sep=""), colnames(data))]) #bandwidth 5nm
  R445_SIRIS <- rowMeans(data[ ,match(paste("rflt_", 444:446, sep=""), colnames(data))]) 
  R470_SIRIS <- rowMeans(data[ ,match(paste("rflt_", 469:471, sep=""), colnames(data))])
  R470to480 <- rowMeans(data[ ,match(paste("rflt_", 470:480, sep=""), colnames(data))])
  R480_ASD <- rowMeans(data[ ,match(paste("rflt_", 476:484, sep=""), colnames(data))])
  R500 <- rowMeans(data[ ,match(paste("rflt_", 499:501, sep=""), colnames(data))])
  R500_SIRIS <- rowMeans(data[ ,match(paste("rflt_", 499:501, sep=""), colnames(data))])
  R500to599 <- rowMeans(data[ ,match(paste("rflt_", 500:599, sep=""), colnames(data))])
  R510 <- rowMeans(data[ ,match(paste("rflt_", 509:511, sep=""), colnames(data))])
  R510_10 <- rowMeans(data[ ,match(paste("rflt_", 506:514, sep=""), colnames(data))])
  R510to520 <- rowMeans(data[ ,match(paste("rflt_", 510:520, sep=""), colnames(data))])
  R513_ASD <- rowMeans(data[ ,match(paste("rflt_", 512:514, sep=""), colnames(data))])
  R520_ASD <- rowMeans(data[ ,match(paste("rflt_", 519:521, sep=""), colnames(data))])
  R520to585 <- rowMeans(data[ ,match(paste("rflt_", 520:585, sep=""), colnames(data))])
  R531_10 <- rowMeans(data[ ,match(paste("rflt_", 527:535, sep=""), colnames(data))])
  R531_6 <- rowMeans(data[ ,match(paste("rflt_", 529:533, sep=""), colnames(data))])
  R534_ASD <- rowMeans(data[ ,match(paste("rflt_", 533:535, sep=""), colnames(data))])
  R540to560 <- rowMeans(data[ ,match(paste("rflt_", 540:560, sep=""), colnames(data))])
  R550 <- rowMeans(data[ ,match(paste("rflt_", 549:551, sep=""), colnames(data))])
  R550_30 <- rowMeans(data[ ,match(paste("rflt_", 535:565, sep=""), colnames(data))])
  R550_ASD <- rowMeans(data[ ,match(paste("rflt_", 549:551, sep=""), colnames(data))])
  R550_CASI <- rowMeans(data[ ,match(paste("rflt_", 547:553, sep=""), colnames(data))])
  R550_HYP <- rowMeans(data[ ,match(paste("rflt_", 543:553, sep=""), colnames(data))]) #B020 548.92 11.02
  R550_SE3 <- rowMeans(data[ ,match(paste("rflt_", 549:551, sep=""), colnames(data))]) #doublons!!
  R550_YARA <- rowMeans(data[ ,match(paste("rflt_", 549:551, sep=""), colnames(data))]) #doublons!!
  R560_10 <- rowMeans(data[ ,match(paste("rflt_", 556:564, sep=""), colnames(data))])
  R560_MSR <- rowMeans(data[ ,match(paste("rflt_", 556:564, sep=""), colnames(data))]) #bandwidth 8.7nm
  R570_10 <- rowMeans(data[ ,match(paste("rflt_", 566:574, sep=""), colnames(data))])
  R570_6 <- rowMeans(data[ ,match(paste("rflt_", 568:572, sep=""), colnames(data))])
  R570_ASD <- rowMeans(data[ ,match(paste("rflt_", 569:571, sep=""), colnames(data))])
  R573_SD <- rowMeans(data[ ,match(paste("rflt_", 567:579, sep=""), colnames(data))]) #bandwidth 14nm
  R584_ASD <- rowMeans(data[ ,match(paste("rflt_", 583:585, sep=""), colnames(data))])
  R600 <- rowMeans(data[ ,match(paste("rflt_", 599:601, sep=""), colnames(data))])
  R600to699 <- rowMeans(data[ ,match(paste("rflt_", 600:699, sep=""), colnames(data))])
  R630 <- rowMeans(data[ ,match(paste("rflt_", 629:631, sep=""), colnames(data))])
  R650_SIRIS <- rowMeans(data[ ,match(paste("rflt_", 649:651, sep=""), colnames(data))])
  R660_MSR <- rowMeans(data[ ,match(paste("rflt_", 656:664, sep=""), colnames(data))]) #bandwidth 9.4nm
  R670 <- rowMeans(data[ ,match(paste("rflt_", 669:671, sep=""), colnames(data))])
  R670_10 <- rowMeans(data[ ,match(paste("rflt_", 666:674, sep=""), colnames(data))])
  R670_5 <- rowMeans(data[ ,match(paste("rflt_", 668:672, sep=""), colnames(data))]) #interpolated ASD
  R670_3 <- rowMeans(data[ ,match(paste("rflt_", 669:671, sep=""), colnames(data))])
  R670_CASI <- rowMeans(data[ ,match(paste("rflt_", 667:673, sep=""), colnames(data))]) #see above Rr_CASI
  R670_SE3 <- rowMeans(data[ ,match(paste("rflt_", 669:671, sep=""), colnames(data))]) #doublons!!
  R670_YARA <- rowMeans(data[ ,match(paste("rflt_", 669:671, sep=""), colnames(data))]) #doublons!!
  R675_SIRIS <- rowMeans(data[ ,match(paste("rflt_", 674:676, sep=""), colnames(data))])
  R678 <- rowMeans(data[ ,match(paste("rflt_", 677:679, sep=""), colnames(data))])
  R680 <- rowMeans(data[ ,match(paste("rflt_", 679:681, sep=""), colnames(data))])
  R680_6 <- rowMeans(data[ ,match(paste("rflt_", 678:682, sep=""), colnames(data))])
  R680_MERIS <- rowMeans(data[ ,match(paste("rflt_", 678:684, sep=""), colnames(data))]) #B08 681.25 7.5
  R680_SIRIS <- rowMeans(data[ ,match(paste("rflt_", 679:681, sep=""), colnames(data))]) #doublons!!
  R680_UNI <- rowMeans(data[ ,match(paste("rflt_", 676:684, sep=""), colnames(data))])
  R681_HYP <- rowMeans(data[ ,match(paste("rflt_", 676:686, sep=""), colnames(data))]) #B033 681.2 10.33
  R683 <- rowMeans(data[ ,match(paste("rflt_", 682:684, sep=""), colnames(data))])
  R690to710 <- rowMeans(data[ ,match(paste("rflt_", 690:710, sep=""), colnames(data))])
  R690to720 <- rowMeans(data[ ,match(paste("rflt_", 690:720, sep=""), colnames(data))])
  R695to740 <- rowMeans(data[ ,match(paste("rflt_", 695:740, sep=""), colnames(data))])
  R698_ASD <- rowMeans(data[ ,match(paste("rflt_", 697:699, sep=""), colnames(data))])
  R700 <- rowMeans(data[ ,match(paste("rflt_", 699:701, sep=""), colnames(data))])
  R700_10 <- rowMeans(data[ ,match(paste("rflt_", 696:704, sep=""), colnames(data))])
  R700_15 <- rowMeans(data[ ,match(paste("rflt_", 693:707, sep=""), colnames(data))])
  R700_ASD <- rowMeans(data[ ,match(paste("rflt_", 699:701, sep=""), colnames(data))]) #doublons!!
  R700_CASI <- rowMeans(data[ ,match(paste("rflt_", 697:703, sep=""), colnames(data))])
  R700_SE3 <- rowMeans(data[ ,match(paste("rflt_", 699:701, sep=""), colnames(data))]) #doublons!!
  R700_YARA <- rowMeans(data[ ,match(paste("rflt_", 699:701, sep=""), colnames(data))]) #doublons!!
  R700_ASD <- rowMeans(data[ ,match(paste("rflt_", 699:701, sep=""), colnames(data))])
  R704_ASD <- rowMeans(data[ ,match(paste("rflt_", 703:705, sep=""), colnames(data))])
  R704 <- rowMeans(data[ ,match(paste("rflt_", 703:705, sep=""), colnames(data))])
  R705_HYP <- rowMeans(data[ ,match(paste("rflt_", 698:706, sep=""), colnames(data))]) #B035 701.55 10.46
  R705 <- rowMeans(data[ ,match(paste("rflt_", 698:706, sep=""), colnames(data))])
  R710 <- rowMeans(data[ ,match(paste("rflt_", 704:706, sep=""), colnames(data))])
  R710_MERIS <- rowMeans(data[ ,match(paste("rflt_", 705:713, sep=""), colnames(data))]) #B09 708.75 10
  R715_ASD <- rowMeans(data[ ,match(paste("rflt_", 714:716, sep=""), colnames(data))])
  R720_10 <- rowMeans(data[ ,match(paste("rflt_", 716:724, sep=""), colnames(data))])
  R720_ASD <- rowMeans(data[ ,match(paste("rflt_", 719:721, sep=""), colnames(data))])
  R720_CASI <- rowMeans(data[ ,match(paste("rflt_", 717:723, sep=""), colnames(data))])
  R724_ASD <- rowMeans(data[ ,match(paste("rflt_", 723:725, sep=""), colnames(data))])
  R726_ASD <- rowMeans(data[ ,match(paste("rflt_", 725:727, sep=""), colnames(data))])
  R730_YARA <- rowMeans(data[ ,match(paste("rflt_", 729:731, sep=""), colnames(data))])
  R734_ASD <- rowMeans(data[ ,match(paste("rflt_", 733:735, sep=""), colnames(data))])
  R737_ASD <- rowMeans(data[ ,match(paste("rflt_", 736:738, sep=""), colnames(data))])
  R740_ASD <- rowMeans(data[ ,match(paste("rflt_", 739:741, sep=""), colnames(data))]) #doublons!!
  R740_YARA <- rowMeans(data[ ,match(paste("rflt_", 739:741, sep=""), colnames(data))])
  R745_ASD <-  rowMeans(data[ ,match(paste("rflt_", 744:746, sep=""), colnames(data))])
  R747_ASD <- rowMeans(data[ ,match(paste("rflt_", 746:748, sep=""), colnames(data))])
  R750 <- rowMeans(data[ ,match(paste("rflt_", 749:751, sep=""), colnames(data))])
  R750_CASI <- rowMeans(data[ ,match(paste("rflt_", 747:753, sep=""), colnames(data))])
  R750_HYP <- rowMeans(data[ ,match(paste("rflt_", 748:756, sep=""), colnames(data))]) #B040 752.43 10.71
  R750_MERIS <- rowMeans(data[ ,match(paste("rflt_", 751:757, sep=""), colnames(data))]) #B10 753.75 7.5
  R750to800 <- rowMeans(data[ ,match(paste("rflt_", 750:800, sep=""), colnames(data))]) #doublons!! see Rn_1
  R760_YARA <- rowMeans(data[ ,match(paste("rflt_", 759:761, sep=""), colnames(data))])
  R760to800 <- rowMeans(data[ ,match(paste("rflt_", 760:800, sep=""), colnames(data))])
  R780_YARA <- rowMeans(data[ ,match(paste("rflt_", 779:781, sep=""), colnames(data))])
  R790_10 <- rowMeans(data[ ,match(paste("rflt_", 786:794, sep=""), colnames(data))])
  R800 <- rowMeans(data[ ,match(paste("rflt_", 799:801, sep=""), colnames(data))])
  R800_10 <- rowMeans(data[ ,match(paste("rflt_", 796:804, sep=""), colnames(data))])
  R800_ASD <- rowMeans(data[ ,match(paste("rflt_", 799:801, sep=""), colnames(data))])
  R800_CASI <- rowMeans(data[ ,match(paste("rflt_", 797:803, sep=""), colnames(data))])
  R800_HYP <- rowMeans(data[ ,match(paste("rflt_", 795:805, sep=""), colnames(data))]) #B045 803.3 11.1
  R800_SIRIS <- rowMeans(data[ ,match(paste("rflt_", 799:801, sep=""), colnames(data))])
  R800_UNI <- rowMeans(data[ ,match(paste("rflt_", 796:804, sep=""), colnames(data))])
  R810_10 <- rowMeans(data[ ,match(paste("rflt_", 806:814, sep=""), colnames(data))])
  R810_MSR <- rowMeans(data[ ,match(paste("rflt_", 805:815, sep=""), colnames(data))]) #bandwidth 11.2nm
  R830 <- rowMeans(data[ ,match(paste("rflt_", 829:831, sep=""), colnames(data))])
  R840_CS <- rowMeans(data[ ,match(paste("rflt_", 770:910, sep=""), colnames(data))])
  R850 <- rowMeans(data[ ,match(paste("rflt_", 849:851, sep=""), colnames(data))])
  R850_5 <- rowMeans(data[ ,match(paste("rflt_", 848:852, sep=""), colnames(data))]) #interpolated ASD
  R860_AVIRIS <- rowMeans(data[ ,match(paste("rflt_", 856:864, sep=""), colnames(data))])
  R874 <- rowMeans(data[ ,match(paste("rflt_", 856:864, sep=""), colnames(data))])
  R880 <- rowMeans(data[ ,match(paste("rflt_", 879:881, sep=""), colnames(data))])
  R893_HYP <- rowMeans(data[ ,match(paste("rflt_", 888:898, sep=""), colnames(data))]) #B075 892.28 11.05
  R900 <- rowMeans(data[ ,match(paste("rflt_", 899:901, sep=""), colnames(data))])
  R900_UNI <- rowMeans(data[ ,match(paste("rflt_", 896:904, sep=""), colnames(data))])
  R900_YARA <- rowMeans(data[ ,match(paste("rflt_", 899:901, sep=""), colnames(data))])
  R920 <- rowMeans(data[ ,match(paste("rflt_", 919:921, sep=""), colnames(data))])
  R920to940 <- rowMeans(data[ ,match(paste("rflt_", 920:940, sep=""), colnames(data))])
  R955 <-  rowMeans(data[ ,match(paste("rflt_", 954:956, sep=""), colnames(data))])
  R960to990 <- rowMeans(data[ ,match(paste("rflt_", 960:990, sep=""), colnames(data))])
  R970 <- rowMeans(data[ ,match(paste("rflt_", 969:971, sep=""), colnames(data))])
  R970_UNI <- rowMeans(data[ ,match(paste("rflt_", 966:974, sep=""), colnames(data))])
  R970_YARA <- rowMeans(data[ ,match(paste("rflt_", 969:971, sep=""), colnames(data))])
  Rb_1 <- rowMeans(data[ ,match(paste("rflt_", 460:480, sep=""), colnames(data))]) #Kaufman and Tanre 1992
  Rb_bb <- rowMeans(data[ ,match(paste("rflt_", 450:520, sep=""), colnames(data))])
  Rb_cb <- rowMeans(data[ ,match(paste("rflt_", 400:520, sep=""), colnames(data))])
  Rb_MODIS3 <- rowMeans(data[ ,match(paste("rflt_", 459:479, sep=""), colnames(data))])
  Re_MERIS <- rowMeans(data[ ,match(paste("rflt_", 700:710, sep=""), colnames(data))])
  Rg_1 <- rowMeans(data[ ,match(paste("rflt_", 530:570, sep=""), colnames(data))]) #MODIS band?
  Rg_bb <- rowMeans(data[ ,match(paste("rflt_", 520:600, sep=""), colnames(data))])
  Rg_cb <- rowMeans(data[ ,match(paste("rflt_", 580:610, sep=""), colnames(data))])
  Rg_MODIS12 <- rowMeans(data[ ,match(paste("rflt_", 546:556, sep=""), colnames(data))])
  Rg_MODIS4 <- rowMeans(data[ ,match(paste("rflt_", 545:565, sep=""), colnames(data))])
  Rg_RE2 <- rowMeans(data[ ,match(paste("rflt_", 520:590, sep=""), colnames(data))]) #RapidEye B2
  Rn_1 <- rowMeans(data[ ,match(paste("rflt_", 845:885, sep=""), colnames(data))]) #Kaufman and Tanre 1992
  Rn_2 <- rowMeans(data[ ,match(paste("rflt_", 750:800, sep=""), colnames(data))]) #doublons!!
  Rn_3 <- rowMeans(data[ ,match(paste("rflt_", 800:900, sep=""), colnames(data))])
  Rn_4 <- rowMeans(data[ ,match(paste("rflt_", 750:900, sep=""), colnames(data))])
  Rn_AVHRR <- rowMeans(data[ ,match(paste("rflt_", 750:1000, sep=""), colnames(data))])
  Rn_bb <- rowMeans(data[ ,match(paste("rflt_", 760:900, sep=""), colnames(data))]) #doublons!!
  Rn_CASI <- rowMeans(data[ ,match(paste("rflt_", 857:863, sep=""), colnames(data))])
  Rn_MERIS <- rowMeans(data[ ,match(paste("rflt_", 750:757, sep=""), colnames(data))])
  Rn_MODIS2 <- rowMeans(data[ ,match(paste("rflt_", 841:876, sep=""), colnames(data))])
  Rn_RE5 <- rowMeans(data[ ,match(paste("rflt_", 760:850, sep=""), colnames(data))]) #RapidEye B5
  Rn_SPOT <- rowMeans(data[ ,match(paste("rflt_", 790:890, sep=""), colnames(data))])
  Rn_TM4 <- rowMeans(data[ ,match(paste("rflt_", 760:900, sep=""), colnames(data))]) #doublons!!
  Rr_1 <- rowMeans(data[ ,match(paste("rflt_", 635:685, sep=""), colnames(data))]) #Kaufman and Tanre 1992 
  Rr_ALI4 <- rowMeans(data[ ,match(paste("rflt_", 630:690, sep=""), colnames(data))]) #doublons!!
  Rr_ASTER2 <- rowMeans(data[ ,match(paste("rflt_", 630:690, sep=""), colnames(data))]) #doublons!!
  Rr_AVHRR <- rowMeans(data[ ,match(paste("rflt_", 580:680, sep=""), colnames(data))])
  Rr_bb <- rowMeans(data[ ,match(paste("rflt_", 630:690, sep=""), colnames(data))]) #doublons!!
  Rr_CASI <- rowMeans(data[ ,match(paste("rflt_", 667:673, sep=""), colnames(data))])
  Rr_cb <- rowMeans(data[ ,match(paste("rflt_", 580:660, sep=""), colnames(data))])
  Rr_MODIS1 <- rowMeans(data[ ,match(paste("rflt_", 620:670, sep=""), colnames(data))])
  Rr_RE3 <- rowMeans(data[ ,match(paste("rflt_", 630:685, sep=""), colnames(data))]) #RapidEye B3
  Rr_SPOT <- rowMeans(data[ ,match(paste("rflt_", 610:690, sep=""), colnames(data))])
  Rr_TM3 <- rowMeans(data[ ,match(paste("rflt_", 630:690, sep=""), colnames(data))]) #doublons!!
  Rre_bb <- rowMeans(data[ ,match(paste("rflt_", 700:730, sep=""), colnames(data))])
  Rre_MERIS <- rowMeans(data[ ,match(paste("rflt_", 704:714, sep=""), colnames(data))])
  Rre_RE4 <- rowMeans(data[ ,match(paste("rflt_", 690:730, sep=""), colnames(data))]) #RapidEye B4
  TM5 <- rowMeans(data[ ,match(paste("rflt_", 1550:1750, sep=""), colnames(data))])
  TM7 <- rowMeans(data[ ,match(paste("rflt_", 2080:2350, sep=""), colnames(data))])
  # A_vis <- rowMeans(data[ ,match(paste("rflt_", 350:700, sep=""), colnames(data))])  #ALBEDO in respective region
  A_NIR <- rowMeans(data[ ,match(paste("rflt_", 700:1000, sep=""), colnames(data))])
  A_SWIR1 <- rowMeans(data[ ,match(paste("rflt_", 1000:1350, sep=""), colnames(data))]) 
  A_SWIR2 <- rowMeans(data[ ,match(paste("rflt_", 1440:1810, sep=""), colnames(data))])
  A_SWIR3 <- rowMeans(data[ ,match(paste("rflt_", 1940:2400, sep=""), colnames(data))])
  R <- rowMeans(data[ ,match(paste("rflt_", 675:685, sep=""), colnames(data))])
  G <- rowMeans(data[ ,match(paste("rflt_", 545:555, sep=""), colnames(data))])
  B <- rowMeans(data[ ,match(paste("rflt_", 495:505, sep=""), colnames(data))])
  
  # INDICES ########################################################################## -
  
  #ARVI
  SI_ARVI <- (Rn_1-(Rr_1-Rb_1+Rr_1))/(Rn_1+(Rr_1-Rb_1+Rr_1)) #with gamma equals 1 here
  
  #EVI/SARVI2
  SI_EVI <- 2.5*(Rn_MODIS2-Rr_MODIS1)/(Rn_MODIS2+6*Rr_MODIS1-7.5*Rb_MODIS3+1)
  
  #NDVI narrow band
  SI_NDVI_nb_ASD <- (R800_ASD-R670_3)/(R800_ASD+R670_3)
  SI_NDVI_nb_CASI <- (R800_CASI-R670_CASI)/(R800_CASI+R670_CASI)
  
  #NDVI MODIS
  SI_NDVI_MODIS <- (Rn_MODIS2-Rr_MODIS1)/(Rn_MODIS2+Rr_MODIS1)
  
  #NDVI broad band
  SI_NDVI_bb <- (Rn_2-Rr_bb)/(Rn_2+Rr_bb)
  SI_NDVI_bb2 <- (Rn_3-Rr_bb)/(Rn_3+Rr_bb)
  SI_NDVI_bb3 <- (Rn_4-Rr_bb)/(Rn_4+Rr_bb)
  
  #SAVI
  L <- 0.5
  SI_SAVI <- (1+L)*(Rn_CASI-Rr_CASI)/(Rn_CASI+Rr_CASI+L)
  
  #NGRDI 
  SI_NGRDI <- (Rg_bb-Rr_bb)/(Rg_bb+Rr_bb)
  
  #WDRVI
  a <- 0.1 #weighting coeff
  SI_WDRVI <- (a*Rn_AVHRR-Rr_AVHRR)/(a*Rn_AVHRR+Rr_AVHRR)
  
  #CAI
  f <- 100  #bands reflectance translation in %
  SI_CAI <-  0.5*(R2020_IRIS*f+R2220_IRIS*f)-R2100_IRIS*f
  
  #RRDI
  SI_RRDI <- (R745_ASD - R740_ASD)/(R740_ASD - R700_ASD)
  
  #mND705
  SI_mND705 <- (R750 - R705)/(R750 + R705 - 2*R445_SIRIS)
  
  #NDTI
  SI_NDTI <- (TM5-TM7)/(TM5+TM7)
  
  #LCA
  SI_LCA <- 100*((ASTER6-ASTER5)/(ASTER6-ASTER8))
  
  #NDRI
  SI_NDRI <- (Rr_TM3-TM7)/(Rr_TM3+TM7)
  
  #ANSI
  SI_ANSI <- (ASTER5-ASTER4)/(ASTER5+ASTER4)
  
  #SINDRI
  SI_SINDRI <- (ASTER6-ASTER7)/(ASTER6+ASTER7)
  
  #NDLI
  SI_NDLI = (log(1/R1754_AVIRIS)-log(1/R1680_AVIRIS))/(log(1/R1754_AVIRIS)+log(1/R1680_AVIRIS))
  
  #SR_dm
  SI_SR_dm <- R780_YARA/R670_YARA
  
  #REIP
  SI_REIP <- 700 + 40*((R670_YARA+R780_YARA)/2-R700_YARA)/(R740_YARA-R700_YARA)
  
  #760/730 ratio
  SI_760_730 <- R760_YARA/R730_YARA
  
  #780/740 ratio
  SI_780_740 <- R780_YARA/R740_YARA
  
  #780/700 ratio  
  SI_780_700	<-R780_YARA/R700_YARA
  
  #780/550 ratio
  SI_780_550 <-R780_YARA/R550_YARA
  
  #970/900 ratio
  SI_970_900 <- R970_YARA/R900_YARA
  
  #MCARI1
  SI_MCARI1 <- 1.2*(2.5*(R800_CASI-Rr_CASI)-1.3*(R800_CASI-R550_CASI))
  
  #MTVI1
  SI_MTVI1 <- 1.2*(1.2*(R800_CASI-R550_CASI)-2.5*(R670_CASI-R550_CASI))
  
  #MCARI2
  SI_MCARI2 <- (1.5*(2.5*(R800_CASI-R670_CASI)-1.3*(R800_CASI-R550_CASI)))/sqrt((2*R800_CASI+1)^2-(6*R800_CASI-5*sqrt(R670_CASI))-0.5)
  
  #MTVI2
  SI_MTVI2 <- (1.5*(1.2*(R800_CASI-R550_CASI)-2.5*(R670_CASI-R550_CASI)))/sqrt((2*R800_CASI+1)^2-(6*R800_CASI-5*sqrt(R670_CASI))-0.5)
  
  #GLI
  SI_GLI <- (2*Rg_cb-Rr_cb-Rb_cb)/(2*Rg_cb+Rr_cb+Rb_cb)
  
  #OSAVI
  SI_OSAVI <- (1+0.16)*(R800_CASI-R670_CASI)/(R800_CASI+R670_CASI+0.16)
  
  #MSAVI
  SI_MSAVI <- 0.5*(2*Rn_SPOT+1-sqrt((2*Rn_SPOT+1)^2-8*(Rn_SPOT-Rr_SPOT)))
  
  #SARVI
  L <- 0.5
  SI_SARVI <- (1+L)*(Rn_1-(Rr_1-Rb_1+Rr_1))/(Rn_1+(Rr_1-Rb_1+Rr_1+L)) #with gamma equals 1 here  #VARIgreen
  SI_VARIgreen <- (Rg_MODIS12-Rr_MODIS1)/(Rg_MODIS12+Rr_MODIS1-Rb_MODIS3)
  
  #VARI700
  SI_VARI700 <- (Re_MERIS-1.7*Rr_MODIS1+0.7*Rb_MODIS3)/(Re_MERIS+2.3*Rr_MODIS1-1.3*Rb_MODIS3)
  
  #VIgreen
  SI_VIgreen <- (Rg_MODIS12-Rr_MODIS1)/(Rg_MODIS12+Rr_MODIS1)
  
  #VI700
  SI_VI700 <- (Re_MERIS-Rr_MODIS1)/(Re_MERIS+Rr_MODIS1)
  
  #SGR
  SI_SGR <- R500to599
  
  #SR_lai
  SI_SR_lai_bb <- Rn_bb/Rr_bb
  SI_SR_lai_cb <- Rn_bb/Rr_cb #no NIR-band for camera
  
  #SLAIDI
  S <- 5 #scaling factor
  SI_SLAIDI1 <- S*(R1050_ASD-R1250_ASD)/(R1050_ASD+R1250_ASD)
  S <- 40 #scaling factor
  SI_SLAIDI2 <- S*R1555_ASD*(R1050_ASD-R1250_ASD)/(R1050_ASD+R1250_ASD)
  
  #SIPI
  SI_SIPI <- (R800_SIRIS-R445_SIRIS)/(R800_SIRIS-R680_SIRIS)
  
  #PSSR
  SI_PSSR1 <- R800_SIRIS/R675_SIRIS
  SI_PSSR2 <- R800_SIRIS/R650_SIRIS
  SI_PSSR3 <- R800_SIRIS/R500_SIRIS
  
  #PSND
  SI_PSND1 <- (R800_SIRIS-R675_SIRIS)/(R800_SIRIS+R675_SIRIS)
  SI_PSND2 <- (R800_SIRIS-R650_SIRIS)/(R800_SIRIS+R650_SIRIS)
  SI_PSND3 <- (R800_SIRIS-R500_SIRIS)/(R800_SIRIS+R500_SIRIS)
  SI_PSND4 <- (R800_SIRIS-R470_SIRIS)/(R800_SIRIS+R470_SIRIS)
  
  SI_PSRI <- (R678-R500)/R750
  
  #NPCI
  SI_NPCI <- (R680_6-R430_6)/(R680_6+R430_6)
  
  #PRI
  SI_PRI <- (R531_10-R570_10)/(R531_10+R570_10)
  
  #PRI570
  SI_PRI570_10 <- (R570_10-R531_10)/(R570_10+R531_10)
  SI_PRI570_6 <- (R570_6-R531_6)/(R570_6+R531_6)
  
  #PRInorm
  SI_PRInorm <- ((R570_10-R531_10)/(R570_10+R531_10))/(((R800_10-R670_10)/sqrt(R800_10+R670_10))*(R700_10/R670_10))
  
  #RGR
  SI_RGR <- R683/R510
  SI_RGR2 <- R600to699/R500to599
  
  #ANTH
  SI_ANTH <- R760to800*(1/R540to560-1/R690to710)
  
  #ARI
  SI_ARI <- 1/(R550*f)-1/(R700*f)
  
  #CARG
  SI_CARG <- R760to800*(1/R510to520-1/R540to560)
  
  #CARRE
  SI_CARRE <- R760to800*(1/R510to520-1/R690to710)
  
  #CRI1
  SI_CRI1 <- 1/R510_10-1/R550_30
  
  #CRI2
  SI_CRI2 <- 1/R510_10-1/R700_15
  
  #CHLG
  SI_CHLG <- R760to800/R540to560
  
  #CHLRE
  SI_CHLRE <- R760to800/R690to720-1
  
  #LCI
  SI_LCI <- (R850-R710)/(R850-R680)
  SI_LCI2 <- (R850-R710)/(R850+R680) #in Pu 2012*
  
  #GNDVI
  SI_GNDVI_HI <- (R750-Rg_1)/(R750+Rg_1) #in Gitelson, Kaufman and Merzlyak 1996
  SI_GNDVI_IRIS <- (R750-R550)/(R750+R550) #in datat 1998
  
  #NDREI
  SI_NDREI <- (R750-R705)/(R750+R705)
  #SI_NDREI2 <- (Rn_bb-Rre_bb)/(Rn_bb+Rre_bb) #in Hunt et al 2011
  
  #NDRE
  SI_NDRE <- (R790_10-R720_10)/(R790_10+R720_10)
  
  #CIG
  SI_CIG <- R750to800/R520to585-1
  #SI_CIG2 <- Rn_bb/Rg_bb-1 #in Hunt et al 2011
  
  #CIRE
  SI_CIRE <- R750to800/R695to740-1
  #SI_CIRE2 <- Rn_bb/Rre_bb-1 #in Hunt et al 2011
  
  #PBI
  SI_PBI <- R810_10/R560_10
  
  #TGI
  SI_TGI_cb <- -0.5*(190*(Rr_cb-Rg_cb)-120*(Rr_cb-Rb_cb))
  SI_TGI_bb <- -0.5*(190*(Rr_bb-Rg_bb)-120*(Rr_bb-Rb_bb))
  SI_TGI_nb <- -0.5*(190*(R670_10-R550_ASD)-120*(R670_10-R480_ASD))
  
  #TCARI
  SI_TCARI <- 3*((R700_CASI-R670_CASI)-0.2*(R700_CASI-R550_CASI)*(R700_CASI/R670_CASI))
  
  #MCARI revised
  SI_MCARI_rev <- ((R750_HYP-R705_HYP)-0.2*(R750_HYP-R550_HYP))/(R750_HYP/R705)
  
  #MSR revised
  SI_MSR_rev <- ((R750_HYP/R705_HYP)-1)/sqrt((R750_HYP/R705_HYP)+1)
  
  #TCARI/OSAVI
  SI_TCARItemp <- 3*((R700_CASI*f-R670_CASI*f)-0.2*(R700_CASI*f-R550_CASI*f)*(R700_CASI*f/(R670_CASI*f)))
  SI_OSAVItemp <- (1+0.16)*(R800_CASI*f-R670_CASI*f)/(R800_CASI*f+R670_CASI*f+0.16)
  SI_TCARI_OSAVI <- SI_TCARItemp/SI_OSAVItemp
  
  #TCARI/OSAVI revised
  TCARItemp <- 3*((R750_HYP-R705_HYP)-0.2*(R750_HYP-R550_HYP)*(R750_HYP/R705_HYP))
  OSAVItemp <- (1+0.16)*(R750_HYP-R705_HYP)/(R750_HYP+R705_HYP+0.16)
  SI_TCARI_OSAVI_rev <- TCARItemp/OSAVItemp
  
  #MCARI/OSAVI revised
  SI_MCARI_OSAVI_rev <- SI_MCARI_rev/OSAVItemp
  
  #MTCI
  SI_MTCI <- (R750_MERIS-R710_MERIS)/(R710_MERIS-R680_MERIS)
  
  #OCAR
  SI_OCAR <- R630/R680
  
  #YCAR
  SI_YCAR <- R600/R680
  
  #MCARI
  SI_MCARI <- ((R700_SE3-R670_SE3)-0.2*(R700_SE3-R550_SE3))*(R700_SE3/R670_SE3)
  
  #TVI
  SI_TVI <- 0.5*(120*(R750_CASI-R550_CASI)-200*(R670_CASI-R550_CASI))
  
  #REM
  SI_REM <- (Rn_MERIS/Rre_MERIS)-1
  #SI_REM2 <- R750/R720-1 #in Chen et al 2010
  
  #GM
  SI_GM <- Rn_MODIS2/Rg_MODIS4-1
  
  #VIopt
  SI_VIopt <- (1+0.45)*((Rn_bb)*2+1)/(Rr_bb+0.45)
  #SI_VIopt2 <- (1+0.45)*((R800)*2+1)/(R670+0.45) #in Chen et al 2010
  
  #RVI1
  SI_RVI1 <- R810_MSR/R660_MSR
  
  #RVI2
  SI_RVI2 <- R810_MSR/R660_MSR
  
  #MCARI/MTVI2
  MCARItemp1 <- (R700-R670-0.2*(R700-R550))*(R700-R670)
  MTVI2temp1 <- (1.5*(1.2*(R800-R550)-2.5*(R670-R550)))/sqrt((2*R800+1)^2-(6*R800-5*sqrt(R670))-0.5)
  MCARItemp2 <- (Rre_RE4-Rr_RE3-0.2*(Rre_RE4-Rg_RE2))*(Rre_RE4-Rr_RE3)
  MTVI2temp2 <- (1.5*(1.2*(Rn_RE5-Rg_RE2)-2.5*(Rr_RE3-Rg_RE2)))/sqrt((2*Rn_RE5+1)^2-(6*R800-5*sqrt(Rr_RE3))-0.5)
  SI_MCARI_MTVI2_ASD <- MCARItemp1/MTVI2temp1
  SI_MCARI_MTVI2_RE <- MCARItemp2/MTVI2temp2
  
  #NDNI
  SI_NDNI <- (log(1/R1510_AVIRIS)-log(1/R1680_AVIRIS))/(log(1/R1510_AVIRIS)+log(1/R1680_AVIRIS))
  
  #DCNI
  n <- 0.03 #soil constant
  SI_DCNI_CASI <- (R720_CASI-R700_CASI)/(R700_CASI-R670_CASI)/(R720_CASI-R670_CASI+n)
  SI_DCNI_ASD <- (R720_ASD-R700_ASD)/(R700_ASD-R670_3)/(R720_ASD-R670_3+n)
  
  #GBNDVI
  SI_GBNDVI <- (R573_SD-R440_SD)/(R573_SD+R440_SD)
  
  #SRWI
  SI_SRWI <- Rn_MODIS2/MODIS5
  
  #WI-WBI
  SI_WI <- R900_UNI/R970_UNI
  
  #WI/NDVI
  SI_WI_NDVI <- (R900_UNI/R970_UNI)/((R800_UNI-R680_UNI)/(R800_UNI+R680_UNI))
  
  #NDWI
  SI_NDWI1 <- (R860_AVIRIS-R1240_AVIRIS)/(R860_AVIRIS+R1240_AVIRIS)
  
  #Ratio at 975
  SI_975 <- (2*R960to990)/(R920to940+R1090to1110)
  
  #Ratio at 1200
  SI_1200 <- (2*R1180to1220)/(R1090to1110+R1265to1285)
  
  #NDMI
  SI_NDMI <- (Rn_TM4-TM5)/(Rn_TM4+TM5)
  
  #LWVI1
  SI_LWVI1 <- (R1094_HYP-R893_HYP)/(R1094_HYP+R893_HYP)
  
  #LWVI2
  SI_LWVI2 <- (R1094_HYP-R1205_HYP)/(R1094_HYP+R1205_HYP)
  
  #LWI
  SI_LWI <- R1300_100/R1450_100
  
  #SIWSI-NDWI1640
  SI_SIWSI <- (Rn_MODIS2-MODIS6)/(Rn_MODIS2+MODIS6)
  
  #NDWI1650
  SI_NDWI1650 <- (R840_CS-R1650_CS)/(R840_CS+R1650_CS)
  
  #NDWI2130
  SI_NDWI2130 <- (Rn_MODIS2-MODIS7)/(Rn_MODIS2+MODIS7)
  
  #PRI
  SI_PRI_10 <- (R531_10-R570_10)/(R531_10+R570_10) #SE590 (10nm)
  SI_PRI_6 <- (R531_6-R570_6)/(R531_6+R570_6) #SE590 (6nm)
  
  #DSWI
  SI_DSWI <- (R800_HYP+R550_HYP)/(R1659_HYP+R681_HYP)
  
  #HI
  SI_HI <- (R534_ASD-R698_ASD)/(R534_ASD+R698_ASD)-R704_ASD/2
  
  #CLSI
  SI_CLSI <- (R698_ASD-R570_ASD)/(R698_ASD+R570_ASD)-R734_ASD
  
  #SBRI
  SI_SBRI <- (R570_ASD-R513_ASD)/(R570_ASD+R513_ASD)-R704_ASD/2
  
  #PMI
  SI_PMI <- (R520_ASD-R584_ASD)/(R520_ASD+R584_ASD)+R724_ASD
  
  #NHI
  SI_NHI_ASD <- (R1100_5-R1200_5)/(R1100_5+R1200_5)
  SI_NHI_ALI <- (ALI6-ALI7)/(ALI6+ALI7)
  
  #Corrected NHI
  SI_CNHI_ASD <- ((R1100_5-R1200_5)/(R1100_5+R1200_5))/((R850_5-R670_5)/(R850_5+R670_5))
  SI_CNHI_ALI <- ((ALI6-ALI7)*(ALI5+Rr_ALI4))/((ALI6+ALI7)*(ALI5-Rr_ALI4))
  
  #FII
  SI_FII <- (R470to480-R360to370)/(R470to480+R360to370)
  
  #NDSVI
  SI_NDSVI <- (TM5-Rr_TM3)/(TM5+Rr_TM3)
  #SI_NDSVI2 <- (ASTER4-ASTER2/ASTER4+ASTER2) #in Pena-Barragan et al 2011
  
  #ALBEDO
  
  # ALBEDO <- (A_vis+A_NIR+A_SWIR1+A_SWIR2+A_SWIR3)/5
  ALBEDO_SWIR <- (A_SWIR1+A_SWIR2+A_SWIR3)/5
  # ALBEDO_VIS <- A_vis*1
  ALBEDO_NIR <- A_NIR*1
  
  ######## SENESCENCE RELEVANT INDICES
  
  #Simple Ratio Index 
  SI_SR <- R800/R670
  
  #Vogelmann Red Edge Index 1
  SI_VOG1 <- R740_ASD/R720_ASD
  
  #Vogelmann Red Edge Index 2
  SI_VOG2 <- (R734_ASD - R747_ASD)/(R715_ASD + R726_ASD)
  
  #Vogelmann Red Edge Index 3
  SI_VOG3 <- (R734_ASD - R747_ASD)/(R715_ASD + R720_ASD)
  
  #R550
  SI_R550 <- R550_ASD
  
  #Gnyli
  SI_GNYLI <- (R900*R1050_ASD-R955*R1220)/(R900*R1050_ASD+R955*R1220)
  
  #NRI
  SI_NRI <- (R874-R1225)/(R874+R1225)
  
  #chormatographic indices (provided by HA)
  #ExG
  SI_Exg_HA <- 2*G-R-B
  
  #NGRDI
  SI_NGRDI_HA <- (G-R)/(G+R)
  
  #GCC
  SI_GCC_HA <- G/(R + G + B)
  
  # CREATE TIBBLE ########################################################################## -
  
  DF <- do.call(cbind.data.frame, mget(ls(pattern = "SI_"))) %>% tibble::as.tibble()
  df_si <- dplyr::bind_cols(data[1:2], DF)
  
  return(df_si)
  
}

#compute NEW spectral indices for wvlt search
f_calc_sen_si <- function(data, wvlt_d, wvlt_n) {
  
  ##wvlt_d: wavelengths to be used as denominators
  ##wvlt_n: wavelengths to be used as numerators
  
  #for each wvlt to be used in the divisor, 
  #extract mean rflt from three adjacent wvlts 
  #and store output in list
  R_denom <- vector("list", length(wvlt_d)) 
  for(i in 1:length(wvlt_d)){
    w <- wvlt_d[i]
    R_denom[[i]] <- rowMeans(data[ ,match(paste("rflt_", (w-1):(w+1), sep=""), colnames(data))])
  }
  
  #name list elements
  names(R_denom) <- paste("R", wvlt_d, sep = "")
  
  #for each wvlt to be used in the zähler, 
  #extract mean rflt from three adjacent wvlts 
  #and store output in list
  R_num <- vector("list", length(wvlt_n)) 
  for(i in 1:length(wvlt_n)){
    w <- wvlt_n[i]
    R_num[[i]] <- rowMeans(data[ ,match(paste("rflt_", (w-1):(w+1), sep=""), colnames(data))])
  }
  
  #name list elements
  names(R_num) <- paste("R", wvlt_n, sep = "")
  
  #extract mean rflt for defined wvlts
  #and add output to list
  R500 <- rowMeans(data[ ,match(paste("rflt_", 499:501, sep=""), colnames(data))])
  R678 <- rowMeans(data[ ,match(paste("rflt_", 677:679, sep=""), colnames(data))])
  R765 <- rowMeans(data[ ,match(paste("rflt_", 764:766, sep=""), colnames(data))])
  
  #------------------- INDICES ----------------------------
  
  #calculate SI
  #PSRI with optimised divisor
  SI0 <- lapply(R_denom, function(x) (R678-R500)/x)
  names(SI0) <- paste("SI0_", names(SI0), sep = "")
  #SR with optimised divisor
  SI1 <- lapply(R_denom, function(x) (R678)/x)
  names(SI1) <- paste("SI1_", names(SI1), sep = "")
  #SR with optimised zahler
  SI2 <- lapply(R_num, function(x) x/(R765))
  names(SI2) <- paste("SI2_", names(SI2), sep = "")
  #PSRI with optimised second band
  SI3 <- lapply(R_num, function(x) (R678-x)/R765)
  names(SI3) <- paste("SI3_", names(SI3), sep = "")
  
  #-------------------  create table ----------------------------
  
  Plot_ID <- data$Plot_ID
  meas_date <- data$meas_date
  
  Result <- do.call("cbind", c(SI0, SI1, SI2, SI3))
  
  d <- data.frame(Plot_ID, meas_date, Result)
  
  return(d)
  
}

#scale spectral indices
f_scale_si <- function(data, sub = NULL) {
  
  #remove measurements taken before heading if required
  if(sub == "post_heading"){
    data <- data[!data$meas_date == "2016-05-26" & !data$meas_date == "2017-05-29", ]
  }
  
  df_sc <- data %>% dplyr::select(-meas_date) %>% nest(-Plot_ID) %>% 
    #scale SI on a plot level
    mutate(data = purrr::map(data, col_mapping)) %>% 
    unnest() %>% 
    dplyr::select(-Plot_ID) %>% 
    #append inverted SI values
    mutate_all(funs(r = revert)) %>% 
    #select original or reversed values
    dplyr::select_if(function(col) col[1] > 5) %>% 
    #complete tibble
    bind_cols(Plot_ID = data$Plot_ID, meas_date = data$meas_date) %>% 
    dplyr::select(Plot_ID, meas_date, everything())
  
  return(df_sc)
  
}  

#bin spectra
f_spc_bin <- function(data, bin_size) {
  
  #Prepare data for binning
  #convert to df and remove info
  data <- data %>% as.data.frame()
  dat <- data %>% dplyr::select(contains("rflt_"))
  #strip rflt from names
  names(dat) <- as.character(gsub("_", "", stringr::str_sub(names(dat), -4, -1)))
  #add rownames
  rownames(dat) <- paste(1:nrow(dat))
  dat <- as.matrix(dat)  
  #Define wavelengths
  wvlt <- as.numeric(colnames(dat))
  
  #Binning
  spc_binned <- prospectr::binning(X = dat, bin.size = bin_size) %>% 
    tibble::as_tibble()
  
  #reassemble tibble
  #correct names
  names(spc_binned) <- paste("rflt_", names(spc_binned), sep = "")
  #add info and convert to tibble
  dat <- data %>% dplyr::select(-contains("rflt_")) %>% 
    dplyr::bind_cols(., spc_binned) %>% tibble::as_tibble()
  
  return(dat)
  
}

#trim spectra
f_spc_trim <- function(data, wb1_lower = 1350, wb1_upper = 1475, 
                       wb2_lower = 1781, wb2_upper = 1990, 
                       final = 2400) {
  
  #Default values similar to Li et al., 2014; Remote Sens.
  
  #bands to drop
  waterband1 <- seq(wb1_lower, wb1_upper, 1)
  waterband2 <- seq(wb2_lower, wb2_upper, 1)
  finalband <- seq(final, as.numeric(stringr::str_sub(names(data[length(data)]), -4, -1)), 1)
  #corresponding variable names
  drop_vars1 <- paste("rflt_", waterband1, sep = "")
  drop_vars2 <- paste("rflt_", waterband2, sep = "")
  drop_vars3 <- paste("rflt_", finalband, sep = "")
  drop_vars <- c(drop_vars1, drop_vars2, drop_vars3)
  
  #clean data
  data <- data[!names(data) %in% drop_vars]
  
  return(data)
  
}

#add sensor information for spectra and vip plotting
add_spc_range <- function(data) {
  
  data <- data %>% 
    mutate(spc_range = ifelse(data$wvlt <= 1350, "visnir", 
                              ifelse(data$wvlt <= 1990, "swir1",
                                     "swir2")))
}

#plot spectra
f_spc_plot <- function(data, by = "meas_date") {
  
  #prepare data
  dat <- data %>%
    tidyr::gather(wvlt, reflectance, contains("rflt_"), factor_key = TRUE) %>% 
    mutate(wvlt = wvlt %>% gsub(x = ., "rflt_", "") %>% as.numeric()) %>% 
    add_spc_range()
  
  #create plot
  Plot <- ggplot(dat) + 
    geom_line(aes(x = wvlt, y = reflectance,  group = interaction(Plot_ID, spc_range), color = Plot_ID), size = 0.2) + 
    scale_y_continuous(limits = c(-0.1, 0.75)) +
    theme_bw() + 
    geom_abline(intercept = 0, slope = 0) +
    facet_wrap( ~ meas_date, ncol=1) 
  
  return(Plot)
  
}

# read senescence rating data
f_sen_read <- function(dir, file_names) {
  
  #list files
  setwd(dir)
  files <- as.list(file_names)
  
  #read files
  d <- lapply(files, function(i){
    x <- read.csv(i)
    x
  })
  
  #add Experiment ID
  #build final tibble
  data <- do.call("rbind", d) %>% mutate(Exp = stringr::str_sub(Plot_ID, 1, 7)) %>% 
    dplyr::select(Exp, everything()) %>% tibble::as_tibble() %>% 
    #transform to dates
    mutate(grading_date = grading_date %>% as.Date(),
           heading_date = heading_date %>% as.Date())
  
  return(data)
  
}

# invert ratings
f_invert_sen <- function(data) {
  
  data <- data %>% mutate(SnsFl0 = SnsFl0 %>% revert(),
                          SnsCnp = SnsCnp %>% revert()) %>% 
    arrange(Plot_ID, grading_date)
  
  return(data)
  
}

# scale ratings
f_scale_sen <- function(data) {
  
  #drop plots which contain missing values in scoring time series
  data <- data %>% group_by(Plot_ID) %>% filter(!any(is.na(SnsCnp))) %>% ungroup() %>% droplevels()
  #select and group
  df_sc <- data %>% arrange(Plot_ID) %>% dplyr::select(Plot_ID, SnsFl0, SnsCnp) %>% 
    nest(-Plot_ID) %>% 
    #scale scorings on a plot level
    mutate(data = purrr::map(data, col_mapping)) %>% 
    unnest() %>% 
    #complete tibble
    bind_cols(data %>% dplyr::select(-SnsCnp, -SnsFl0, -Plot_ID), .) %>% 
    #reorder columns
    dplyr::select(Plot_ID, everything())
  
  return(df_sc)
  
}

# read *.asd files
f_spc_read <- function(dir, Exp) {
  
  #files have to be named measDate_"ASD"_PlotID_"Can"_rep.asd
  
  #load data
  setwd(dir)
  fnames <- list.files(dir, pattern = ".asd")
  ASD <- prospectr::readASD(fnames, 'binary', 'list')
  rflt <- sapply(ASD, "[[", "reflectance")
  wvlt <- sapply(ASD, "[[", "wavelength")
  names <- sapply(fnames, basename)
  
  #transpose to wide format
  df <- cbind(wvlt[,1], rflt)
  df <- t(df)[-1,] %>% as.data.frame()
  
  #rename columns
  colnames(df) <- unique(paste("rflt", wvlt, sep = "_"))
  
  #extract info from spectra name
  pn0 <- names %>% strsplit(., "_", fixed = TRUE) %>% unlist()
  Plot_ID <- grep(Exp, pn0, value = TRUE)
  replicate <- stringr::str_sub(names,-5,-5)
  meas_date <- stringr::str_sub(strptime(grep("201", names, value = TRUE), "%Y%m%d"), 1, 10) %>% as.Date()
  
  #create data frame
  df <- data.frame(Plot_ID, meas_date, replicate, df)
  
  #drop reference from data frame
  df <- df[!grepl("Ref", df$Plot_ID),]
  df$Plot_ID <- df$Plot_ID %>% droplevels()
  df <- tibble::as_tibble(df)
  
  return(df)
  
}

# perform continuum removal
f_cont_rem <- function(data, interpol = "linear", method = "substraction") {
  
  #Prepare required format
  #convert to df and remove info
  data <- data %>% as.data.frame() #to enable rownames
  dat <- data %>% dplyr::select(contains("rflt_"))
  #strip rflt from names
  names(dat) <- as.character(gsub("_", "", str_sub(names(dat), -4, -1)))
  #add rownames
  rownames(dat) <- paste(1:nrow(dat))
  dat <- as.matrix(dat)
  # Define wavelengths
  wvlt <- as.numeric(colnames(dat))
  
  # continuum removal
  spc_cr <- prospectr::continuumRemoval(X = dat, wvlt, type = "R", 
                                        interpol = interpol, method = method) %>% tibble::as_tibble()
  
  #reassemble tibble
  #correct names
  names(spc_cr) <- paste("rflt_", names(spc_cr), sep = "")
  #add info and convert to tibble
  dat <- data %>% dplyr::select(-contains("rflt_")) %>% 
    dplyr::bind_cols(spc_cr) %>% tibble::as_tibble()
  
  return(dat)
  
}

# define matching dates and join data sets
is.date <- function(x) inherits(x, 'Date')

f_match_join <- function(spc, sen, gddah, matches) {
  
  #check whether dates are correctly specified (as.Date)!
  if(!is.date(sen$grading_date)){sen$grading_date = as.Date(sen$grading_date)}
  if(!is.date(gddah$heading_date)){gddah$heading_date = as.Date(gddah$heading_date)}
  
  #matching dates lookup table
  match_dates <- matches %>% transmute(scor_date = scor %>% as.Date(),
                                       meas_date = meas %>%  as.Date(),
                                       ref_date = reference_date %>% as.Date())
  
  #assign scoring dates their measurement date, 
  #if possible (max 1 day difference)
  sen_data <- sen %>% full_join(., match_dates, by = c("grading_date" = "scor_date"))
  
  #convert gddah data to tibble
  gddah_data <- gddah %>% tibble::as_tibble() %>% 
    mutate(meas_date = meas_date %>% as.Date(),
           heading_date = heading_date %>% as.Date())
  
  #add gddah data to spectral data
  spc_data <- right_join(gddah, spc, by = c("Plot_ID", "meas_date")) %>% 
    #Plots for which heading date is missing must be excluded  
    filter(!is.na(heading_date)) 
  
  #join spectral and scoring data
  data <- full_join(spc_data, sen_data) %>% 
    #reorder columns
    dplyr::select(-contains("rflt_"), everything()) %>% 
    dplyr::select(Exp, Plot_ID, everything()) %>% 
    arrange(Plot_ID, ref_date)
  
}

# extract errors and senescence dynamics parameters
get_errors_and_dynpars <- function(data, method){
  
  if (method != "lin"){
    
    if (method == "log"){
      
      #fit logistic model to scoring and SVI values
      m1 <- nls(formula = as.formula("scoring ~ A + C/(1+exp(-b*(grading_GDDAH-M)))"), data = data, 
                start = list(A = 10, C = -10, b = 0.01, M = 675), na.action = na.exclude)
      m2 <- nls(formula = as.formula("value ~ A + C/(1+exp(-b*(meas_GDDAH-M)))"), data = data, 
                start = list(A = 10, C = -10, b = 0.01, M = 675), na.action = na.exclude)
      
    } else if (method == "cgom"){
      
      #fit constrained Gompertz model to scoring and SVI values
      m1 <- nls(formula = as.formula("scoring ~10*exp(-exp(-b*(grading_GDDAH-M)))"), data = data, 
                start = list(b = 0.01, M = 675), na.action = na.exclude)
      m2 <- nls(formula = as.formula("value ~ 10*exp(-exp(-b*(meas_GDDAH-M)))"), data = data, 
                start = list(b = 0.01, M = 675), na.action = na.exclude)
      
    } else if (method == "gom"){
      
      #fit flexible Gompertz model to scoring and SVI values
      m1 <- nls(formula = as.formula("scoring ~ A+C*exp(-exp(-b*(grading_GDDAH-M)))"), data = data, 
                start = list(A = 0, C = 11, b = 0.01, M = 675), na.action = na.exclude)
      m2 <- nls(formula = as.formula("value ~ A+C*exp(-exp(-b*(meas_GDDAH-M)))"), data = data, 
                start = list(A = 0, C = 11, b = 0.01, M = 675), na.action = na.exclude)
    }
    
    #get min and max and create sequence of values to predict
    r1 <- range(data$grading_GDDAH, na.rm = TRUE)
    xNew1 <- seq(r1[1],r1[2],length.out = 1000)
    
    #create predictions
    y <- predict(m1, list(grading_GDDAH = xNew1))
    y2 <- predict(m2, list(meas_GDDAH = xNew1))
    
    #create output identical to linear interpolation (see below)
    l <- list(xNew1, y)
    ll <- list(xNew1, y2)
    names(l) <- names(ll) <- c("x", "y")
    
  } else if (method == "lin"){
    
    x <- data$meas_GDDAH
    y <- data$value
    
    #linear interpolation of SVI values
    l <- approx(x = x, y = y, 
                xout = seq(round(min(x, na.rm = TRUE),0), 
                           round(max(x, na.rm = TRUE), 0)))
    
    x <- data$grading_GDDAH
    y <- data$scoring
    
    #linear interpolation of Scorings
    ll <- approx(x = x, y = y, 
                 xout = seq(round(min(x, na.rm = TRUE), 0), 
                            round(max(x, na.rm = TRUE), 0), 1))
  }
  
  #Calculate error
  d1 <- as.data.frame(do.call("cbind", l))
  d2 <- as.data.frame(do.call("cbind", ll))
  d3 <- merge(d1, d2, by = "x")
  d3$y <- abs(d3$y.x - d3$y.y)
  d3$y.x[is.na(d3$y.x)] <- 0
  d3$y.y[is.na(d3$y.y)] <- 0
  f1 <- approxfun(d3$x, d3$y.x-d3$y.y)     # piecewise linear function
  f2 <- function(x) abs(f1(x))
  Error <- integrate(f2, min(d3$x), max(d3$x), subdivisions = 2000)$value
  
  #Extract senescence dynamics parameters
  onsen_SI <- l[[1]][which(l[[2]] < 8)[1]]
  onsen_Trait <- ll[[1]][which(ll[[2]] < 8)[1]]
  midsen_SI <- l[[1]][which(l[[2]] < 5)[1]]
  midsen_Trait <- ll[[1]][which(ll[[2]] < 5)[1]]
  endsen_SI <- l[[1]][which(l[[2]] < 2)[1]]
  endsen_Trait <- ll[[1]][which(ll[[2]] < 2)[1]]
  tsen_SI <- endsen_SI - onsen_SI
  tsen_Trait <- endsen_Trait - onsen_Trait
  
  #create data frame
  func_out <- do.call(rbind, Map(tibble, "onsen_SI" = onsen_SI, "midsen_SI" = midsen_SI,
                                 "endsen_SI" = endsen_SI, "tsen_SI" = tsen_SI,
                                 "onsen_Trait" = onsen_Trait, "midsen_Trait" = midsen_Trait,
                                 "endsen_Trait" = endsen_Trait, "tsen_Trait" = tsen_Trait)) %>% 
    #add bias
    mutate(d_onsen = onsen_SI - onsen_Trait,
           d_midsen = midsen_SI - midsen_Trait,
           d_endsen = endsen_SI - endsen_Trait,
           d_tsen = tsen_SI - tsen_Trait) %>% 
    #add error
    mutate(Error = Error)
  
  return(func_out)
  
}

# extract senescence dynamics parameters from spectral indices
get_dynpars_SI <- function(data, method){
  
  x <- data$meas_GDDAH
  y <- data$value
  
  #linear interpolation of SVI values
  l <- approx(x = x, y = y, 
              xout = seq(round(min(x, na.rm = TRUE),0), 
                         round(max(x, na.rm = TRUE), 0)))
  
  #Extract senescence dynamics parameters
  onsen_SI <- l[[1]][which(l[[2]] < 8)[1]]
  midsen_SI <- l[[1]][which(l[[2]] < 5)[1]]
  endsen_SI <- l[[1]][which(l[[2]] < 2)[1]]
  tsen_SI <- endsen_SI - onsen_SI

  #create data frame
  func_out <- do.call(rbind, Map(tibble, "onsen_SI" = onsen_SI, "midsen_SI" = midsen_SI,
                                 "endsen_SI" = endsen_SI, "tsen_SI" = tsen_SI)) 
 
  return(func_out)
  
}


#============================================================================================ -

# HELPER FUNCTIONS ----

#calculate correlation and p-values
do_cor_test <- function(data, x, y, use = "pairwise.complete.obs", 
                        method = "pearson", return = "estimate") {
  cor <- cor.test(data %>% pull(x), data %>% pull(y), 
                  use = use, method = method) %>% 
    broom::tidy(.) %>% pull(return)
  return(cor)
}
  

# Contrained Gompertz equation
Gompertz_constrained <- function(b, M, tempsum) {
  grenness_decay <- 10*exp(-exp(-b*(tempsum-M)))
  return(grenness_decay)
}

# Logistic equation
logistic <- function(A, C, b, M, tempsum) {
  grenness_decay <- A + C/(1+exp(-b*(tempsum-M)))
  return(grenness_decay)
}

# Flexible Gompertz equation
Gompertz_flex <- function(A, C, b, M, tempsum) {
  grenness_decay <- A + C*exp(-exp(-b*(grading_GDDAH-M)))
  return(grenness_decay)
}

# Linear interpolation
lin_approx <- function(data){
  data <- as.data.frame(data)
  out <- approx(data[,"grading_GDDAH"], data[,"Score"], xout = seq(round(min(data[,"grading_GDDAH"], na.rm = TRUE), 0),
                                                                   round(max(data[,"grading_GDDAH"], na.rm = TRUE), 0), 1))
  names(out) <- c("grading_GDDAH", ".fitted")
  out <- list(out)
  return(out)
}

# Extraction of parameters from fits
extract_pars <- function(data){
  onsen <- data[which(data[2] < 8)[1], "grading_GDDAH"]
  midsen <- data[which(data[2] < 5)[1], "grading_GDDAH"]
  endsen <- data[which(data[2] < 2)[1], "grading_GDDAH"]
  tsen <- endsen - onsen
  pars <- cbind(onsen, midsen, endsen, tsen)
  names(pars) <- c("onsen", "midsen", "endsen", "tsen")
  return(pars)
}

## column mapping for scaling
col_mapping <- function(d) {
  purrr::map_df(d, function(X) (X - min(X, na.rm = TRUE))/diff(range(X, na.rm = TRUE))*10)
}

## revert si values
revert <- function(x){10 - x}