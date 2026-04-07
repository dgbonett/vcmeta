# replicate.mean2 returns valid matrix

    Code
      res
    Output
                            Estimate        SE      t     p        LL       UL     df
      Original:                 5.80 0.7889312  7.352 0.000  4.228624 7.371376  75.75
      Follow-up:                6.10 0.6346075  9.612 0.000  4.845913 7.354087 147.65
      Original - Follow-up:    -0.30 1.0124916 -0.296 0.767 -1.974571 1.374571 169.16
      Average:                  5.95 0.5062458 11.753 0.000  4.950627 6.949373 169.16

# replicate.mean.ps returns valid matrix

    Code
      res
    Output
                            Estimate       SE     t     p        LL       UL   df
      Original:                15.29 2.154344 7.097 0.000 10.780906 19.79909 19.0
      Follow-up:                7.57 1.460664 5.183 0.000  4.659564 10.48044 74.0
      Original - Follow-up:     7.72 2.602832 2.966 0.005  3.332885 12.10712 38.4
      Average:                 11.43 1.301416 8.783 0.000  8.796322 14.06368 38.4

# replicate.stdmean2 returns valid matrix

    Code
      res
    Output
                            Estimate     SE      LL     UL
      Original:               1.6280 0.2595  1.1353 2.1524
      Follow-up:              1.5617 0.2044  1.1690 1.9704
      Original - Follow-up:   0.0742 0.3303 -0.4691 0.6176
      Average:                1.5949 0.1652  1.2711 1.9186

# replicate.stdmean.ps returns valid matrix

    Code
      res
    Output
                            Estimate     SE     LL     UL
      Orginal:                1.0890 0.2292 0.6697 1.5680
      Follow-up:              0.4605 0.0959 0.2757 0.6516
      Original - Follow-up:   0.6552 0.2484 0.2466 1.0638
      Average:                0.7748 0.1242 0.5313 1.0182

# replicate.cor returns valid matrix

    Code
      res
    Output
                            Estimate      SE     z     p     LL     UL
      Original:                0.598 0.07321 6.589 0.000 0.4355 0.7228
      Follow-up:               0.324 0.06377 4.819 0.000 0.1940 0.4428
      Original - Follow-up:    0.274 0.09709 2.633 0.008 0.1065 0.4265
      Average:                 0.461 0.04854 7.635 0.000 0.3725 0.5412

# replicate.cor.gen returns valid matrix

    Code
      res
    Output
                            Estimate      SE     z     p      LL     UL
      Original:                0.454 0.17000 2.287 0.022  0.0699 0.7209
      Follow-up:               0.318 0.09800 3.022 0.003  0.1152 0.4953
      Original - Follow-up:    0.136 0.19622 0.667 0.505 -0.2154 0.4237
      Average:                 0.386 0.09811 3.409 0.001  0.1961 0.5480

# replicate.gen returns valid matrix

    Code
      res
    Output
                            Estimate        SE     z     p         LL        UL
      Original:                0.782 0.2100000 3.724 0.000  0.3704076 1.1935924
      Follow-up:               0.650 0.1540000 4.221 0.000  0.3481655 0.9518345
      Original - Follow-up:    0.132 0.2604151 0.507 0.612 -0.2963446 0.5603446
      Average:                 0.716 0.1302075 5.499 0.000  0.4607979 0.9712021

---

    Code
      res
    Output
                            Estimate        SE     z     p         LL        UL
      Original:                0.782 0.2100000 3.724 0.000  0.3704076 1.1935924
      Follow-up:               0.650 0.1540000 4.221 0.000  0.3481655 0.9518345
      Original - Follow-up:    0.132 0.2604151 0.507 0.612 -0.2963446 0.5603446
      Average:                 0.716 0.1302075 5.499 0.000  0.4607979 0.9712021

# replicate.oddsratio returns valid matrix

    Code
      res
    Output
                              Estimate         SE     z     p          LL        UL
      Original:             0.11904762 0.10805233 1.102 0.271 -0.09273105 0.3308263
      Follow-up:            0.09677419 0.07965047 1.215 0.224 -0.05933787 0.2528863
      Original - Follow-up: 0.02359056 0.13542107 0.174 0.862 -0.19915727 0.2463384
      Average:              0.11015594 0.06771053 1.627 0.104 -0.02255427 0.2428661

# replicate.prop2 returns valid matrix

    Code
      res
    Output
                             Estimate        SE      z     p exp(Estimate)   exp(LL)
      Original:            1.39000000 0.3020000  4.603 0.000     4.0148501 2.2212961
      Follow-up:           1.48000000 0.2060000  7.184 0.000     4.3929457 2.9336501
      Original/Follow-up: -0.06273834 0.3655681 -0.172 0.864     0.9391892 0.5147653
      Average:             0.36067292 0.1827840  1.973 0.048     1.4342943 1.0024257
                           exp(UL)
      Original:           7.256583
      Follow-up:          6.578144
      Original/Follow-up: 1.713551
      Average:            2.052222

# replicate.slope returns valid matrix

    Code
      res
    Output
                            Estimate       SE     t     p        LL       UL    df
      Original:                23.40 5.160000 4.535 0.000 13.007227 33.79277  45.0
      Follow-up:               18.50 4.480000 4.129 0.000  9.592560 27.40744  85.0
      Original - Follow-up:     4.90 6.833447 0.717 0.475 -6.438743 16.23874 106.4
      Average:                 20.95 3.416724 6.132 0.000 14.176310 27.72369 106.4

# replicate.spear returns valid matrix

    Code
      res
    Output
                            Estimate      SE     z     p     LL     UL
      Original:                0.598 0.07948 5.315 0.000 0.4199 0.7318
      Follow-up:               0.324 0.06542 4.571 0.000 0.1905 0.4457
      Original - Follow-up:    0.274 0.10294 3.438 0.001 0.0948 0.4342
      Average:                 0.461 0.05147 9.968 0.000 0.3670 0.5457

# replicate.mean1 returns valid matrix

    Code
      res
    Output
                            Estimate        SE        LL        UL     df
      Original:                21.90 0.6039950 20.678305 23.121695 39.000
      Follow-up:               25.20 0.4595708 24.284285 26.115715 74.000
      Original - Follow-up:    -3.30 0.7589567 -4.562527 -2.037473 82.630
      Average:                 23.55 0.3794784 22.795183 24.304817 82.633

# replicate.prop1 returns valid matrix

    Code
      res
    Output
                               Estimate         SE          LL         UL
      Original:              0.07565789 0.01516725  0.04593064 0.10538515
      Follow-up:             0.09158416 0.01435033  0.06345803 0.11971029
      Original - Follow-up: -0.01670456 0.02065098 -0.05067239 0.01726328
      Average:               0.08119996 0.01032549  0.06096237 0.10143755

# replicate.propratio2 returns valid matrix

    Code
      res
    Output
                            Estimate        SE      z     p exp(Estimate)   exp(LL)
      Original:            0.2682640 0.2463597  1.089 0.276     1.3076923 0.8068705
      Follow-up:           0.3735135 0.3082711  1.212 0.226     1.4528302 0.7939881
      Original/Follow-up: -0.1052495 0.3946190 -0.267 0.789     0.9000999 0.4703209
      Average:             0.3208887 0.1973095  1.626 0.104     1.3783522 0.9362893
                           exp(UL)
      Original:           2.119373
      Follow-up:          2.658372
      Original/Follow-up: 1.722611
      Average:            2.029132

# replicate.prop.ps returns valid matrix

    Code
      res
    Output
                               Estimate         SE     z     p          LL         UL
      Original:             0.106557377 0.03440159 3.097 0.002  0.03913151 0.17398325
      Follow-up:            0.103174603 0.02358274 4.375 0.000  0.05695329 0.14939592
      Original - Follow-up: 0.003852359 0.04097037 0.094 0.925 -0.06353791 0.07124263
      Average:              0.105511837 0.02048519 5.151 0.000  0.06536161 0.14566206

# replicate.agree example

    Code
      res
    Output
                            Estimate      SE      LL     UL
      Original:                 0.70 0.07252  0.5309 0.8152
      Follow-up:                0.60 0.05662  0.4773 0.6992
      Original - Follow-up:     0.10 0.09160 -0.0584 0.2429
      Average:                  0.65 0.04580  0.5504 0.7299

# replicate.cronbach example

    Code
      res
    Output
                            Estimate      SE      LL     UL
      Original:                0.883 0.01831  0.8436 0.9152
      Follow-up:               0.869 0.01442  0.8387 0.8952
      Original - Follow-up:    0.014 0.02331 -0.0334 0.0582
      Average:                 0.876 0.01172  0.8519 0.8977

# replicate.meanratio2 example

    Code
      res
    Output
                              Estimate         SE       t       p     df
      Original:             0.30766736 0.04188600  7.3454 0.00000  76.65
      Follow-up:            0.27715566 0.02928438  9.4643 0.00000 140.91
      Original - Follow-up: 0.03051171 0.05110785  0.5970 0.55141 150.35
      Average:              0.29241151 0.02555393 11.4429 0.00000 150.35
                            exp(Estimate)   exp(LL)  exp(UL)
      Original:                  1.360248 1.2513908 1.478576
      Follow-up:                 1.319372 1.2451576 1.398009
      Original - Follow-up:      1.030982 0.9473616 1.121983
      Average:                   1.339654 1.2736927 1.409032

# replicate.meanratio.ps example

    Code
      res
    Output
                             Estimate         SE      t       p    df exp(Estimate)
      Original:             0.1952087 0.02655112 7.3522 0.00000 19.00      1.215565
      Follow-up:            0.0934960 0.01839400 5.0830 0.00000 74.00      1.098006
      Original - Follow-up: 0.1017127 0.03230018 3.1490 0.00313 39.29      1.107065
      Average:              0.1443523 0.01615009 8.9382 0.00000 39.29      1.155291
                             exp(LL)  exp(UL)
      Original:             1.149856 1.285028
      Follow-up:            1.058492 1.138996
      Original - Follow-up: 1.048437 1.168972
      Average:              1.118170 1.193645

# replicate.pbcor example

    Code
      res
    Output
                            Estimate      SE      LL      UL
      Original:               0.4753 0.05847  0.3487  0.5781
      Follow-up:              0.6016 0.03365  0.5292  0.6617
      Original - Follow-up:  -0.1262 0.06746 -0.2664 -0.0005
      Average:                0.5384 0.03373  0.4691  0.6012

