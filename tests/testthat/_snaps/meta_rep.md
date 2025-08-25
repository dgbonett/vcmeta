# replicate.mean2 returns valid matrix

    Code
      res
    Output
                            Estimate        SE          t            p        LL
      Original:                 5.80 0.7889312  7.3517180 1.927969e-10  4.228624
      Follow-up:                6.10 0.6346075  9.6122408 0.000000e+00  4.845913
      Original - Follow-up:    -0.30 1.0124916 -0.2962988 7.673654e-01 -1.974571
      Average:                  5.95 0.5062458 11.7531843 0.000000e+00  4.950627
                                  UL        df
      Original:             7.371376  75.75255
      Follow-up:            7.354087 147.64728
      Original - Follow-up: 1.374571 169.16137
      Average:              6.949373 169.16137

# replicate.mean.ps returns valid matrix

    Code
      res
    Output
                            Estimate       SE        t            p        LL
      Original:                15.29 2.154344 7.097288 9.457592e-07 10.780906
      Follow-up:                7.57 1.460664 5.182575 1.831197e-06  4.659564
      Original - Follow-up:     7.72 2.602832 2.966000 5.166213e-03  3.332885
      Average:                 11.43 1.301416 8.782740 1.010232e-10  8.796322
                                  UL       df
      Original:             19.79909 19.00000
      Follow-up:            10.48044 74.00000
      Original - Follow-up: 12.10712 38.40002
      Average:              14.06368 38.40002

# replicate.stdmean2 returns valid matrix

    Code
      res
    Output
                              Estimate        SE        LL        UL
      Original:             1.62803662 0.2594668  1.135349 2.1524396
      Follow-up:            1.56170447 0.2044452  1.168967 1.9703776
      Original - Follow-up: 0.07422178 0.3303344 -0.469130 0.6175736
      Average:              1.59487055 0.1651672  1.271149 1.9185923

# replicate.stdmean.ps returns valid matrix

    Code
      res
    Output
                             Estimate         SE        LL        UL
      Orginal:              1.0890300 0.22915553 0.6697353 1.5680085
      Follow-up:            0.4604958 0.09590506 0.2756687 0.6516096
      Original - Follow-up: 0.6552328 0.24841505 0.2466264 1.0638392
      Average:              0.7747629 0.12420752 0.5313206 1.0182052

# replicate.cor returns valid matrix

    Code
      res
    Output
                            Estimate         SE        z            p        LL
      Original:                0.598 0.07320782 6.589418 4.708045e-09 0.4355043
      Follow-up:               0.324 0.06376782 4.819037 2.865955e-06 0.1939787
      Original - Follow-up:    0.274 0.09708614 2.633335 8.455096e-03 0.1065496
      Average:                 0.461 0.04854307 7.634998 2.264855e-14 0.3725367
                                   UL
      Original:             0.7227538
      Follow-up:            0.4428347
      Original - Follow-up: 0.4265016
      Average:              0.5411607

# replicate.cor.gen returns valid matrix

    Code
      res
    Output
                            Estimate         SE         z            p          LL
      Original:                0.454 0.17000000 2.2869806 0.0221969560  0.06991214
      Follow-up:               0.318 0.09800000 3.0215123 0.0025151541  0.11522137
      Original - Follow-up:    0.136 0.19622436 0.6671281 0.5046902807 -0.21543667
      Average:                 0.386 0.09811218 3.4089419 0.0006521538  0.19606750
                                   UL
      Original:             0.7208577
      Follow-up:            0.4953353
      Original - Follow-up: 0.4237240
      Average:              0.5480170

# replicate.gen returns valid matrix

    Code
      res
    Output
                            Estimate        SE         z            p         LL
      Original:                0.782 0.2100000 3.7238095 1.962390e-04  0.3704076
      Follow-up:               0.650 0.1540000 4.2207792 2.434593e-05  0.3481655
      Original - Follow-up:    0.132 0.2604151 0.5068831 6.122368e-01 -0.2963446
      Average:                 0.716 0.1302075 5.4989141 3.821373e-08  0.4607979
                                   UL
      Original:             1.1935924
      Follow-up:            0.9518345
      Original - Follow-up: 0.5603446
      Average:              0.9712021

---

    Code
      res
    Output
                            Estimate        SE         z            p         LL
      Original:                0.782 0.2100000 3.7238095 1.962390e-04  0.3704076
      Follow-up:               0.650 0.1540000 4.2207792 2.434593e-05  0.3481655
      Original - Follow-up:    0.132 0.2604151 0.5068831 6.122368e-01 -0.2963446
      Average:                 0.716 0.1302075 5.4989141 3.821373e-08  0.4607979
                                   UL
      Original:             1.1935924
      Follow-up:            0.9518345
      Original - Follow-up: 0.5603446
      Average:              0.9712021

# replicate.oddsratio returns valid matrix

    Code
      res
    Output
                              Estimate         SE         z         p          LL
      Original:             0.11904762 0.10805233 1.1017590 0.2705665 -0.09273105
      Follow-up:            0.09677419 0.07965047 1.2149858 0.2243715 -0.05933787
      Original - Follow-up: 0.02359056 0.13542107 0.1742016 0.8617070 -0.19915727
      Average:              0.11015594 0.06771053 1.6268656 0.1037656 -0.02255427
                                   UL
      Original:             0.3308263
      Follow-up:            0.2528863
      Original - Follow-up: 0.2463384
      Average:              0.2428661

# replicate.prop2 returns valid matrix

    Code
      res
    Output
                             Estimate        SE          z            p   exp(LL)
      Original:            1.39000000 0.3020000  4.6026490 4.171509e-06 2.2212961
      Follow-up:           1.48000000 0.2060000  7.1844660 6.747936e-13 2.9336501
      Original/Follow-up: -0.06273834 0.3655681 -0.1716188 8.637372e-01 0.5147653
      Average:             0.36067292 0.1827840  1.9732190 4.847061e-02 1.0024257
                           exp(UL)
      Original:           7.256583
      Follow-up:          6.578144
      Original/Follow-up: 1.713551
      Average:            2.052222

# replicate.slope returns valid matrix

    Code
      res
    Output
                            Estimate       SE         t            p        LL
      Original:                23.40 5.160000 4.5348837 4.250869e-05 13.007227
      Follow-up:               18.50 4.480000 4.1294643 8.465891e-05  9.592560
      Original - Follow-up:     4.90 6.833447 0.7170612 4.749075e-01 -6.438743
      Average:                 20.95 3.416724 6.1316052 1.504129e-08 14.176310
                                  UL       df
      Original:             33.79277  45.0000
      Follow-up:            27.40744  85.0000
      Original - Follow-up: 16.23874 106.4035
      Average:              27.72369 106.4035

# replicate.spear returns valid matrix

    Code
      res
    Output
                            Estimate         SE        z            p         LL
      Original:                0.598 0.07948367 5.315140 1.065752e-07 0.41985966
      Follow-up:               0.324 0.06541994 4.570582 4.863705e-06 0.19049455
      Original - Follow-up:    0.274 0.10294378 3.437975 5.860809e-04 0.09481418
      Average:                 0.461 0.05147189 9.967944 0.000000e+00 0.36695230
                                   UL
      Original:             0.7317733
      Follow-up:            0.4457384
      Original - Follow-up: 0.4342171
      Average:              0.5457190

# replicate.mean1 returns valid matrix

    Code
      res
    Output
                            Estimate        SE        LL        UL       df
      Original:                21.90 0.6039950 20.678305 23.121695 39.00000
      Follow-up:               25.20 0.4595708 24.284285 26.115715 74.00000
      Original - Follow-up:    -3.30 0.7589567 -4.562527 -2.037473 82.63282
      Average:                 23.55 0.3794784 22.795183 24.304817 82.63282

# replicate.prop1 returns valid matrix

    Code
      res
    Output
                               Estimate         SE          LL         UL
      Original:              0.07565789 0.01516725  0.04593064 0.10538515
      Follow-up:             0.09158416 0.01435033  0.06345803 0.11971029
      Original - Follow-up: -0.01670456 0.02065098 -0.05067239 0.01726328
      Average:               0.08119996 0.01032549  0.06096237 0.10143755

# replicate.prop.ps returns valid matrix

    Code
      res
    Output
                               Estimate         SE          z            p
      Original:             0.106557377 0.03440159 3.09745539 1.951898e-03
      Follow-up:            0.103174603 0.02358274 4.37500562 1.214294e-05
      Original - Follow-up: 0.003852359 0.04097037 0.09402793 9.250870e-01
      Average:              0.105511837 0.02048519 5.15064083 2.595979e-07
                                     LL         UL
      Original:              0.03913151 0.17398325
      Follow-up:             0.05695329 0.14939592
      Original - Follow-up: -0.06353791 0.07124263
      Average:               0.06536161 0.14566206

# replicate.agree example

    Code
      res
    Output
                            Estimate         SE          LL        UL
      Original:                 0.70 0.07252105  0.53093828 0.8152156
      Follow-up:                0.60 0.05661961  0.47726289 0.6992077
      Original - Follow-up:     0.10 0.09159681 -0.05844824 0.2428784
      Average:                  0.65 0.04579840  0.55040374 0.7299302

# replicate.cronbach example

    Code
      res
    Output
                            Estimate         SE          LL         UL
      Original:                0.883 0.01830958  0.84356871 0.91522517
      Follow-up:               0.869 0.01442263  0.83874629 0.89523760
      Original - Follow-up:    0.014 0.02330779 -0.03336284 0.05820123
      Average:                 0.876 0.01172239  0.85187755 0.89774525

