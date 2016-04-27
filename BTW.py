#### BTW-FF atom types and properties #####
BTW_atoms = {
#FF_num  at_num  at_mass   valance vdW_rad[A]  epsilo[kcal/mol] H-bond  charge   atom_type  description
"21" :(  1      , 1.008   , 1.0    , 1.62     ,   0.02          , 0.923 ,  0.622   ), #  H         H-Oi                            
"75" :(  8      , 15.9994 , 4.0    , 1.82     ,   0.059         , 0     , -1.242   ), #  O         O-H
"170":(  8      , 15.9994 , 2.0    , 1.82     ,   0.059         , 0     , -1.0908  ), #  O         O-Carboxylate                                
"171":(  8      , 15.9994 , 4.0    , 1.82     ,   0.059         , 0     , -1.1145  ), #  O         O-inorganic
"172":(  30     , 65.38   , 4.0    , 2.29     ,   0.276         , 0     ,  1.281   ), #  Zn        Zn
"192":(  40     , 91.224  , 8.0    , 3.52     ,   0.367         , 0     ,  2.601   ), #  Zr        Zr
"185":(  29     , 63.546  , 5.0    , 2.29     ,   0.276         , 0     ,  1.0358  ), #  Cu        Cu                               
"902":(  6      , 12.0    , 3.0    , 1.96     ,   0.056         , 0     , -0.0114  ), #  C         Calpha
"903":(  6      , 12      , 3.0    , 1.96     ,   0.056         , 0     , -0.0124  ), #  C         C-doublephenolligand
"913":(  6      , 12.0    , 3.0    , 1.94     ,   0.056         , 0     ,  1.5398  ), #  C         Cacid
"912":(  6      , 12.0    , 3.0    , 1.96     ,   0.056         , 0     , -0.0228  ), #  C         Cbenzene
"915":(  1      , 1.008   , 1.0    , 1.62     ,   0.02          , 0.923 ,  0.1582  )  #  H         Hbenzene        
}
####   BONDs in BTW-FF ####
BTW_bonds = {
#FF_type    k[mdyne]     r[A]
"21_75"  :(   3.630      ,   0.989),
"75_192" :(   5.500      ,   2.276),
"170_172":(   3.665      ,   2.009),
"170_192":(   5.821      ,   2.338),
"170_185":(   5.091      ,   1.969),
"170_913":(   5.999      ,   1.299), 
"171_172":(   4.329      ,   2.039),
"171_192":(   5.809      ,   2.192),
"185_185":(   4.349      ,   2.422),
"902_912":(   4.500      ,   1.389), 
"902_913":(   5.299      ,   1.485), 
"903_903":(   5.899      ,   1.465),
"903_912":(   5.999      ,   1.389),
"912_912":(   4.500      ,   1.389), 
"912_915":(   5.150      ,   1.101) 
}
#### ANGLES in BTW-FF ####
BTW_angles = {
# at1_atcen_at2 k[mdyne/rad^2], Theta1[degree] , Theta2[degree], Theta3[degree], Ksb1 , Ksb2 , Kss
"21_75_192"  :(     2.099,            116.848  ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 ),
"75_192_75"  :(     2.099,            123.230  ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 ),
"75_192_170" :(     2.099,             89.658  ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 ),
"75_192_171" :(     2.099,             71.110  ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 ),
"170_172_170":(     1.000,            110.103  ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 ),   
"170_172_171":(     3.000,            113.584  ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 ),   
"170_192_170":(     2.099,             73.103  ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 ),
"170_192_171":(     2.099,             84.318  ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 ),
"170_185_185":(     4.299,             87.822  ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 ),
#"170_185_170":(     0.05,            175.945   ,    0.0        ,    0.0        , 0.0  , 0.0  , 0.0 ),# Fourier equation used instead
"170_913_170":(     2.867,            126.299  ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 ), 
"170_913_902":(     1.867,            117.082  ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 ),     
"171_192_171":(     2.099,             91.479  ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 ), 
"172_170_913":(     3.022,            130.606  ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 ),   
"172_171_172":(     1.198,            110.992  ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 ), 
"185_170_913":(     3.322,            120.962  ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 ),
"192_75_192" :(     2.099,            103.406  ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 ),
"192_170_913":(     2.099,            139.820  ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 ),
"192_171_192":(     2.099,            118.408  ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 ),
"902_912_902":(      0.06,            121.797  ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 ), # Added from MM3 
"902_912_912":(     0.060,            121.582  ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 ),
"902_912_915":(     0.090,            119.859  ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 ), 
"902_913_170":(     1.867,            117.082  ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 ),     
"903_903_912":(     5.00 ,            122.690  ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 ),
"903_912_912":(     5.00 ,            122.904  ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 ),
"903_912_915":(     5.00 ,            120.00   ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 ),
"912_902_912":(     0.000,            119.406  ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 ),
"912_902_913":(     0.360,            121.797  ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 ),
"912_903_912":(     5.00 ,            117.621  ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 ),
"912_912_915":(     0.49 ,            120.0    ,    0.0        ,     0.0       , 0.0  ,  0.0 , 0.0 )  #  Added from MM3 
}
#### DIHEDRALs in BTW-FF
BTW_dihedrals = {
#at1_at2_at3_at4     k1  , t1  , n1 ,  k2    ,  t2    , n2 , k3  , t3  , n3 , k4  , t4  , n4
"170_913_902_912":(  0.0 , 0.0 ,  1 ,  2.5   ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ), 
"915_912_902_913":(  0.0 , 0.0 ,  1 ,  1.999 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"913_902_912_912":(  0.0 , 0.0 ,  1 ,  8.030 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"912_902_912_912":(  0.0 , 0.0 ,  1 ,  8.030 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"912_902_912_915":(  0.0 , 0.0 ,  1 ,  8.030 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"902_912_912_915":(  0.0 , 0.0 ,  1 ,  8.030 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"902_912_912_902":(  0.0 , 0.0 ,  1 ,  8.030 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"903_903_912_915":(  0.0 , 0.0 ,  1 ,  6.999 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"903_903_912_912":(  0.0 , 0.0 ,  1 ,  6.9   ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"912_903_903_912":(  0.0 , 0.0 ,  1 ,  6.9   ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ), 
"912_903_912_912":(  0.0 , 0.0 ,  1 ,  4.930 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"912_903_912_915":(  0.0 , 0.0 ,  1 ,  4.930 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"915_912_912_903":(  0.0 , 0.0 ,  1 ,  4.930 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"902_912_912_903":(  0.0 , 0.0 ,  1 ,  5.930 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"902_912_912_912":(  0.0 , 0.0 ,  1 ,  8.030 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"171_172_170_913":(  0.0 , 0.0 ,  1 ,  4.690 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ), 
"172_170_913_170":(  0.0 , 0.0 ,  1 ,  2.176 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ), 
"170_172_170_913":(  0.0 , 0.0 ,  1 ,  0.860 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ), 
"172_170_913_902":(  0.0 , 0.0 ,  1 ,  0.072 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"170_172_171_172":(  0.0 , 0.0 ,  1 ,  1.000 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"170_185_170_913":(  0.0 , 0.0 ,  1 ,  0.860 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"185_170_913_902":(  0.0 , 0.0 ,  1 ,  0.072 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"170_913_170_185":(  0.0 , 0.0 ,  1 ,  5.805 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"913_170_185_185":(  0.0 , 0.0 ,  1 ,  0.850 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"170_185_185_170":(  0.0 , 0.0 ,  1 ,  2.071 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"170_185_185_170":(  0.0 , 0.0 ,  1 ,  2.071 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"171_192_170_913":(  0.0 , 0.0 ,  1 ,  2.064 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"192_170_913_170":(  0.0 , 0.0 ,  1 ,  2.017 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"915_912_902_913":(  0.0 , 0.0 ,  1 ,  1.999 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"170_192_170_913":(  0.0 , 0.0 ,  1 ,  0.860 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"192_170_913_902":(  0.0 , 0.0 ,  1 ,  0.072 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"170_192_171_192":(  0.0 , 0.0 ,  1 ,  1.000 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"913_170_192_75" :(  0.0 , 0.0 ,  1 ,  5.000 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"21_75_192_170"  :(  0.0 , 0.0 ,  1 ,  5.000 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"21_75_191_170"  :(  0.0 , 0.0 ,  1 ,  5.000 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"192_75_192_170" :(  0.0 , 0.0 ,  1 ,  5.000 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"192_75_192_75"  :(  0.0 , 0.0 ,  1 ,  5.000 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"192_75_192_171" :(  0.0 , 0.0 ,  1 ,  5.000 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"21_75_192_75"   :(  0.0 , 0.0 ,  1 ,  5.000 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"192_75_192_171" :(  0.0 , 0.0 ,  1 ,  5.000 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"192_171_192_75" :(  0.0 , 0.0 ,  1 ,  5.000 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"192_171_192_171":(  0.0 , 0.0 ,  1 ,  5.000 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"21_75_192_171"  :(  0.0 , 0.0 ,  1 ,  5.000 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"192_75_192_171" :(  0.0 , 0.0 ,  1 ,  5.000 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"192_171_192_75" :(  0.0 , 0.0 ,  1 ,  5.000 ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ),
"915_912_912_915":(  0.0 , 0.0 ,  1 ,  11.5  ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ), # Added from MM3
"902_912_902_913":(  0.0 , 0.0 ,  1 ,  8.03  ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ), # Added from MM3
"902_912_902_912":(  0.0 , 0.0 ,  1 ,  8.03  ,  180.0 ,  2 , 0.0 , 0.0 , 3  , 0.0 , 0.0 , 4 ), # Added from MM3
}
#### OUT-OF-PLANE bending in BTW-FF
BTW_opbends = {
#at1_at2_at3_at4   K_opb, phi , Ka1  ,  ka2  , ka3 
###    """
###          H
###         /
###    C = C 
###         \
###          C
###    """
"902_912_912_915":( 0.0 , 0.0 , 0.24 , 0.300 , 0.0 ),#  BTW-FF coefficient is 0.0 while in MM3 is 0.2 
"902_912_915_912":( 0.0 , 0.0 , 0.24 , 0.300 , 0.0 ),#  BTW-FF coefficient is 0.0 while in MM3 is 0.2
"902_912_902_915":( 0.0 , 0.0 , 0.24 , 0.300 , 0.0 ),#  BTW-FF coefficient is 0.0 while in MM3 is 0.2
"902_912_915_902":( 0.0 , 0.0 , 0.24 , 0.300 , 0.0 ),#  BTW-FF coefficient is 0.0 while in MM3 is 0.2
"903_912_915_912":( 0.2 , 0.0 , 0.24 , 0.300 , 0.0 ),#  Added from MM3
"903_912_912_915":( 0.2 , 0.0 , 0.24 , 0.300 , 0.0 ),#  Added from MM3
"912_912_915_903":( 0.2 , 0.0 , 0.24 , 0.300 , 0.0 ),#  Added from MM3
"912_912_903_915":( 0.2 , 0.0 , 0.24 , 0.300 , 0.0 ),#  Added from MM3
"912_912_902_915":( 0.2 , 0.0 , 0.24 , 0.300 , 0.0 ),#  Added from MM3
"912_912_915_902":( 0.2 , 0.0 , 0.24 , 0.300 , 0.0 ),#  Added from MM3
"915_912_912_902":( 0.11, 0.0 , 0.24 , 0.300 , 0.0 ),#  Added from MM3
"915_912_902_912":( 0.11, 0.0 , 0.24 , 0.300 , 0.0 ),#  Added from MM3
"915_912_902_902":( 0.11, 0.0 , 0.24 , 0.300 , 0.0 ),#  Added from MM3 
"915_912_912_903":( 0.11, 0.0 , 0.24 , 0.300 , 0.0 ),#  Added from MM3 
"915_912_903_912":( 0.11, 0.0 , 0.24 , 0.300 , 0.0 ),#  Added from MM3 
###   """
###         O
###        /
###   C = C
###        \
###         O
###   """
"170_913_170_902":( 1.5 , 0.0 , 0.00 , 0.0   , 0.0 ),#  
"170_913_902_170":( 1.5 , 0.0 , 0.00 , 0.0   , 0.0 ),#  
"902_913_170_170":( 0.0 , 0.0 , 0.00 , 0.0   , 0.0 ),# BTW-FF coefficient is 0.0 while in MM3 is 0.2  
###   ----------- ########   -----------
### O |       C | ######## C |       C |
###  \|      /  | ##    ##  \|      /  |
###   |C = C    | ## Or ##   |C = C    |
###  /|      \  | ##    ##  /|      \  |
### O |       C | ######## C |       C |
###   ----------- ########   -----------
"912_902_912_913":( 0.2 , 0.0 , 0.24 , 0.300 , 0.0 ),# Added from MM3
"912_902_913_912":( 0.2 , 0.0 , 0.24 , 0.300 , 0.0 ),# Added from MM3
"913_902_912_912":( 0.2 , 0.0 , 0.24 , 0.300 , 0.0 ),# Added from MM3
"903_903_912_912":( 0.2 , 0.0 , 0.24 , 0.300 , 0.0 ),# Added from MM3
"912_903_912_903":( 0.2 , 0.0 , 0.24 , 0.300 , 0.0 ),# Added from MM3
"912_903_903_912":( 0.2 , 0.0 , 0.24 , 0.300 , 0.0 ),# Added from MM3 
#
# UIO special opbend
#
"192_171_192_192":( 2.0 , 0.0 , 0.00 , 0.00  , 0.0 ) 
}


BTW_charges = {
"Cu Paddlewheel_185"    :(  1.0358 ),  #  Cu         Cu                               
"Cu Paddlewheel_170"    :( -1.0908 ),  #  O          O-Carboxylate                                
"Cu Paddlewheel_902"    :( -0.0114 ),  #  C          Calpha
"Cu Paddlewheel_913"    :(  1.5398 ),  #  C          Cacid
"Cu Paddlewheel_912"    :( -0.0228 ),  #  C          Cbenzene
"Cu Paddlewheel_915"    :(  0.1582 ),  #  H          Hbenzene     
"Zn4O_170"              :( -1.1513 ),  #  O       O-Carboxylate                                
"Zn4O_171"              :( -1.1145 ),  #  O       O-inorganic
"Zn4O_902"              :( -0.0081 ),  #  C       Calpha
"Zn4O_913"              :(  1.4972 ),  #  C       Cacid
"Zn4O_912"              :( -0.0536 ),  #  C       Cbenzene
"Zn4O_915"              :(  0.1259 ),  #  H       Hbenzene                        
"Zn4O_172"              :(  1.281  ),  #  Zn      Zn
"IRMOF10_170"           :( -1.1630 ),  #  O       O-Carboxylate                                
"IRMOF10_902"           :( -0.0279 ),  #  C       Calpha
"IRMOF10_913"           :(  1.5377 ),  #  C       Cacid
"IRMOF10_912"           :( -0.0460 ),  #  C       Cbenzene
"IRMOF10_915"           :(  0.1047 ),  #  H       Hbenzene                        
"IRMOF10_172"           :(  1.2954 ),  #  Zn      Zn
"IRMOF10_171"           :( -1.2144 ),  #  O       O-inorganic
"IRMOF10_903"           :( -0.0124 ),  #  C       C-doublephenolligand
"Zr_UiO_170"            :( -1.181  ),  #  O       O-Carboxylate                                
"Zr_UiO_902"            :( -0.056  ),  #  C       Calpha
"Zr_UiO_913"            :(  1.576  ),  #  C       Cacid
"Zr_UiO_912"            :( -0.058  ),  #  C       Cbenzene
"Zr_UiO_915"            :(  0.129  ),  #  H       Hbenzene                        
"Zr_UiO_75"             :( -1.242  ),  #  O       O-H
"Zr_UiO_192"            :(  2.601  ),  #  Zr      Zr
"Zr_UiO_21"             :(  0.622  ),  #  H       H-Oi                            
"Zr_UiO_171"            :( -1.189  ),  #  O       O-inorganic
"Zr_UiO_903"            :( -0.035  ),  #  C       C-doublephenolligand     
"TFF_171"               :( -1.186  ),  #  O       O-inorganic             
"TFF_172"               :(  1.291  ),  #  Zn      Zn     
"TFF_170"               :( -1.154  ),  #  O       O-Carboxylate          
"TFF_913"               :(  1.539  ),  #  C       Cacid        
"TFF_912"               :( -0.050  ),  #  C       Cbenzene   
"TFF_902"               :( -0.008  ),  #  C       Calpha
"TFF_915"               :(  0.118  ),  #  H       Hbenzene
"TFF_192"               :(  2.605  ),  #  Zr      Zr    
"TFF_75"                :( -1.243  ),  #  O       O-H
"TFF_21"                :(  0.622  ),  #  H       O-H
"TFF_903"               :( -0.0124 ),  #  C       C-doublephenolligand     !!!!!! temporary!
}

