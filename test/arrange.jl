using Plasm

V = [0.7650325 0.3990061 1.0036906 0.6376642 0.9702703 0.6042439 1.2089284 0.842902 0.303951 -0.0226949 0.5933759 0.26673 0.8571261 0.5304803 1.146551 0.8199052 -0.2683782 0.2400275 -0.1672349 0.3411708 -0.2478461 0.2605596 -0.1467028 0.3617028 0.707404 0.2149479 0.929584 0.4371279 0.8197257 0.3272696 1.0419057 0.5494496 0.4157712 0.315036 0.3160769 0.2153417 1.1038847 1.0031495 1.0041904 0.9034552 -0.109785 0.0036575 -0.0037144 0.1097281 0.1699003 0.2833428 0.2759709 0.3894134 0.2240651 -0.252741 0.4878986 0.0110925 0.3993298 -0.0774763 0.6631632 0.1863571 0.476173 -0.0281847 0.2858746 -0.2184831 0.5537047 0.049347 0.3634063 -0.1409514 0.3339188 0.3990749 0.3981511 0.4633072 0.7112856 0.7764417 0.7755179 0.840674; -0.2546268 -0.0159687 0.1610265 0.3996846 -0.3123348 -0.0736767 0.1033185 0.3419766 0.1349662 0.4243911 0.758342 1.0477668 -0.0202849 0.26914 0.6030909 0.8925158 0.9984893 1.0996326 0.5307374 0.6318806 0.7982245 0.8993677 0.3304725 0.4316158 -0.1315806 0.0905994 0.3729568 0.5951368 -0.1554784 0.0667016 0.349059 0.571239 0.4426001 0.3429058 1.1327847 1.0330904 0.5279999 0.4283056 1.2181845 1.1184902 0.2388182 0.3448888 0.5042387 0.6103093 0.0951344 0.201205 0.3605549 0.4666255 0.5018578 0.7656913 1.0079403 1.2717738 0.4577868 0.7216202 0.9638693 1.2277028 0.6651987 0.4749004 1.175287 0.9849886 0.6792641 0.4889658 1.1893524 0.999054 0.0840742 0.1483065 0.4596064 0.5238387 0.0090638 0.0732961 0.3845961 0.4488283; 0.5727381 0.7779758 0.5150301 0.7202678 1.0058694 1.2111072 0.9481615 1.1533992 0.4631022 1.0162774 0.3078512 0.8610263 0.8709765 1.4241516 0.7157254 1.2689005 0.4944629 0.514995 0.2941981 0.3147301 0.9725845 0.9931166 0.7723197 0.7928517 0.050891 0.1632127 0.0269932 0.1393149 0.5906185 0.7029402 0.5667207 0.6790424 0.3752028 1.0633164 0.4606026 1.1487162 0.4883108 1.1764243 0.5737106 1.2618241 0.0147364 0.2944217 -0.1289474 0.1507379 -0.0442141 0.2354712 -0.1878979 0.0917874 0.1918858 0.3671505 0.1478148 0.3230794 0.7350339 0.9102986 0.6909629 0.8662276 -0.1563218 -0.0787901 -0.1422564 -0.0647247 0.3825588 0.4600906 0.3966242 0.4741559 0.4042854 0.7816522 0.3292751 0.7066419 0.351897 0.7292638 0.2768866 0.6542534]

CV = [[1, 2, 3, 4, 5, 6, 7, 8], [9, 10, 11, 12, 13, 14, 15, 16], [17, 18, 19, 20, 21, 22, 23, 24], [25, 26, 27, 28, 29, 30, 31, 32], [33, 34, 35, 36, 37, 38, 39, 40], [41, 42, 43, 44, 45, 46, 47, 48], [49, 50, 51, 52, 53, 54, 55, 56], [57, 58, 59, 60, 61, 62, 63, 64], [65, 66, 67, 68, 69, 70, 71, 72]]
FV = [[1, 2, 3, 4], [5, 6, 7, 8], [1, 2, 5, 6], [3, 4, 7, 8], [1, 3, 5, 7], [2, 4, 6, 8], [9, 10, 11, 12], [13, 14, 15, 16], [9, 10, 13, 14], [11, 12, 15, 16], [9, 11, 13, 15], [10, 12, 14, 16], [17, 18, 19, 20], [21, 22, 23, 24], [17, 18, 21, 22], [19, 20, 23, 24], [17, 19, 21, 23], [18, 20, 22, 24], [25, 26, 27, 28], [29, 30, 31, 32], [25, 26, 29, 30], [27, 28, 31, 32], [25, 27, 29, 31], [26, 28, 30, 32], [33, 34, 35, 36], [37, 38, 39, 40], [33, 34, 37, 38], [35, 36, 39, 40], [33, 35, 37, 39], [34, 36, 38, 40], [41, 42, 43, 44], [45, 46, 47, 48], [41, 42, 45, 46], [43, 44, 47, 48], [41, 43, 45, 47], [42, 44, 46, 48], [49, 50, 51, 52], [53, 54, 55, 56], [49, 50, 53, 54], [51, 52, 55, 56], [49, 51, 53, 55], [50, 52, 54, 56], [57, 58, 59, 60], [61, 62, 63, 64], [57, 58, 61, 62], [59, 60, 63, 64], [57, 59, 61, 63], [58, 60, 62, 64], [65, 66, 67, 68], [69, 70, 71, 72], [65, 66, 69, 70], [67, 68, 71, 72], [65, 67, 69, 71], [66, 68, 70, 72]]
EV = [[1, 2], [3, 4], [5, 6], [7, 8], [1, 3], [2, 4], [5, 7], [6, 8], [1, 5], [2, 6], [3, 7], [4, 8], [9, 10], [11, 12], [13, 14], [15, 16], [9, 11], [10, 12], [13, 15], [14, 16], [9, 13], [10, 14], [11, 15], [12, 16], [17, 18], [19, 20], [21, 22], [23, 24], [17, 19], [18, 20], [21, 23], [22, 24], [17, 21], [18, 22], [19, 23], [20, 24], [25, 26], [27, 28], [29, 30], [31, 32], [25, 27], [26, 28], [29, 31], [30, 32], [25, 29], [26, 30], [27, 31], [28, 32], [33, 34], [35, 36], [37, 38], [39, 40], [33, 35], [34, 36], [37, 39], [38, 40], [33, 37], [34, 38], [35, 39], [36, 40], [41, 42], [43, 44], [45, 46], [47, 48], [41, 43], [42, 44], [45, 47], [46, 48], [41, 45], [42, 46], [43, 47], [44, 48], [49, 50], [51, 52], [53, 54], [55, 56], [49, 51], [50, 52], [53, 55], [54, 56], [49, 53], [50, 54], [51, 55], [52, 56], [57, 58], [59, 60], [61, 62], [63, 64], [57, 59], [58, 60], [61, 63], [62, 64], [57, 61], [58, 62], [59, 63], [60, 64], [65, 66], [67, 68], [69, 70], [71, 72], [65, 67], [66, 68], [69, 71], [70, 72], [65, 69], [66, 70], [67, 71], [68, 72]]

V,CVs,FVs,EVs = testarrangement(V,FV,EV)
show_exploded(V,CVs,FVs,EVs)