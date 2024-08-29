using Plasm
include("../book/exploder.jl")

n = 10
V = [0.9101582 0.6379717 0.6592827 0.3870961 1.5062844 1.2340979 1.2554089 0.9832224 0.0843469 -0.0881351 0.2383533 0.0658713 0.3897303 0.2172483 0.5437367 0.3712547 -0.036538 -0.1151385 0.0390101 -0.0395905 0.2265768 0.1479763 0.3021248 0.2235243 0.7490819 0.6994406 0.6997791 0.6501378 1.1691151 1.1194738 1.1198123 1.070171 0.072093 0.3745675 -0.0689181 0.2335564 0.1464 0.4488745 0.0053889 0.3078634 0.1641349 0.1334626 0.1944844 0.1638121 0.3716781 0.3410058 0.4020276 0.3713553 0.0803545 0.6958805 0.147548 0.7630741 0.0877337 0.7032598 0.1549272 0.7704533 -0.1387411 -0.4444038 0.134104 -0.1715587 0.4015458 0.0958831 0.6743909 0.3687282 -0.1846788 -0.1108464 -0.2580307 -0.1841983 0.4551466 0.528979 0.3817947 0.4556271 0.304116 -0.1708791 0.0209969 -0.4539983 0.5142812 0.0392861 0.2311621 -0.2438331; 0.4028344 0.1519588 1.039914 0.7890385 0.5563971 0.3055215 1.1934767 0.9426012 0.022095 0.1761015 0.3624505 0.516457 -0.0625644 0.091442 0.277791 0.4317975 -0.0292294 0.0463186 0.2398719 0.3154199 -0.0839278 -0.0083798 0.1851735 0.2607215 0.6842707 0.6349679 1.1049786 1.0556758 0.7278259 0.6785231 1.1485337 1.0992309 0.5895819 0.4485708 0.4271526 0.2861415 0.8553432 0.7143321 0.6929139 0.5519028 0.7177194 0.748069 0.9259049 0.9562544 0.6917613 0.7221108 0.8999467 0.9302962 0.6810943 0.7482878 0.0802826 0.1474761 0.5471088 0.6143023 -0.0537029 0.0134906 0.5269135 0.7997587 1.1293121 1.4021573 0.3770618 0.6499069 0.9794604 1.2523055 0.6535834 0.5802315 1.292451 1.2190991 0.7352899 0.661938 1.3741575 1.3008056 0.2619841 -0.0211351 0.7783923 0.4952731 0.3177727 0.0346535 0.8341809 0.5510618; -0.1870511 0.4090752 -0.0334884 0.5626379 0.1497612 0.7458875 0.3033239 0.8994502 -0.0859219 0.2194615 -0.1705814 0.134802 0.1292543 0.4346377 0.0445948 0.3499782 0.4366695 0.6997843 0.3819711 0.6450859 0.5309756 0.7940904 0.4762772 0.739392 0.3690667 0.7890998 0.4126218 0.832655 0.4238204 0.8438536 0.4673756 0.8874087 0.6111982 0.6855052 0.8769595 0.9512665 0.8130543 0.8873613 1.0788156 1.1531226 0.3601404 0.5676836 0.3341822 0.5417255 0.3946086 0.6021518 0.3686504 0.5761937 0.1449814 0.1523606 0.0109959 0.0183751 0.7494938 0.756873 0.6155083 0.6228875 -0.07578 0.4645069 -0.2256318 0.3146551 0.3055579 0.8458448 0.1557061 0.695993 0.5562169 1.1960423 0.6379235 1.2777488 0.4917517 1.1315771 0.5734582 1.2132836 -0.2326003 -0.0224351 -0.1768117 0.0333535 0.3175492 0.5277144 0.3733378 0.583503];
CV = [[1, 2, 3, 4, 5, 6, 7, 8], [9, 10, 11, 12, 13, 14, 15, 16], [17, 18, 19, 20, 21, 22, 23, 24], [25, 26, 27, 28, 29, 30, 31, 32], [33, 34, 35, 36, 37, 38, 39, 40], [41, 42, 43, 44, 45, 46, 47, 48], [49, 50, 51, 52, 53, 54, 55, 56], [57, 58, 59, 60, 61, 62, 63, 64], [65, 66, 67, 68, 69, 70, 71, 72], [73, 74, 75, 76, 77, 78, 79, 80]];
FV = [[1, 2, 3, 4], [5, 6, 7, 8], [1, 2, 5, 6], [3, 4, 7, 8], [1, 3, 5, 7], [2, 4, 6, 8], [9, 10, 11, 12], [13, 14, 15, 16], [9, 10, 13, 14], [11, 12, 15, 16], [9, 11, 13, 15], [10, 12, 14, 16], [17, 18, 19, 20], [21, 22, 23, 24], [17, 18, 21, 22], [19, 20, 23, 24], [17, 19, 21, 23], [18, 20, 22, 24], [25, 26, 27, 28], [29, 30, 31, 32], [25, 26, 29, 30], [27, 28, 31, 32], [25, 27, 29, 31], [26, 28, 30, 32], [33, 34, 35, 36], [37, 38, 39, 40], [33, 34, 37, 38], [35, 36, 39, 40], [33, 35, 37, 39], [34, 36, 38, 40], [41, 42, 43, 44], [45, 46, 47, 48], [41, 42, 45, 46], [43, 44, 47, 48], [41, 43, 45, 47], [42, 44, 46, 48], [49, 50, 51, 52], [53, 54, 55, 56], [49, 50, 53, 54], [51, 52, 55, 56], [49, 51, 53, 55], [50, 52, 54, 56], [57, 58, 59, 60], [61, 62, 63, 64], [57, 58, 61, 62], [59, 60, 63, 64], [57, 59, 61, 63], [58, 60, 62, 64], [65, 66, 67, 68], [69, 70, 71, 72], [65, 66, 69, 70], [67, 68, 71, 72], [65, 67, 69, 71], [66, 68, 70, 72], [73, 74, 75, 76], [77, 78, 79, 80], [73, 74, 77, 78], [75, 76, 79, 80], [73, 75, 77, 79], [74, 76, 78, 80]];
EV = [[1, 2], [3, 4], [5, 6], [7, 8], [1, 3], [2, 4], [5, 7], [6, 8], [1, 5], [2, 6], [3, 7], [4, 8], [9, 10], [11, 12], [13, 14], [15, 16], [9, 11], [10, 12], [13, 15], [14, 16], [9, 13], [10, 14], [11, 15], [12, 16], [17, 18], [19, 20], [21, 22], [23, 24], [17, 19], [18, 20], [21, 23], [22, 24], [17, 21], [18, 22], [19, 23], [20, 24], [25, 26], [27, 28], [29, 30], [31, 32], [25, 27], [26, 28], [29, 31], [30, 32], [25, 29], [26, 30], [27, 31], [28, 32], [33, 34], [35, 36], [37, 38], [39, 40], [33, 35], [34, 36], [37, 39], [38, 40], [33, 37], [34, 38], [35, 39], [36, 40], [41, 42], [43, 44], [45, 46], [47, 48], [41, 43], [42, 44], [45, 47], [46, 48], [41, 45], [42, 46], [43, 47], [44, 48], [49, 50], [51, 52], [53, 54], [55, 56], [49, 51], [50, 52], [53, 55], [54, 56], [49, 53], [50, 54], [51, 55], [52, 56], [57, 58], [59, 60], [61, 62], [63, 64], [57, 59], [58, 60], [61, 63], [62, 64], [57, 61], [58, 62], [59, 63], [60, 64], [65, 66], [67, 68], [69, 70], [71, 72], [65, 67], [66, 68], [69, 71], [70, 72], [65, 69], [66, 70], [67, 71], [68, 72], [73, 74], [75, 76], [77, 78], [79, 80], [73, 75], [74, 76], [77, 79], [78, 80], [73, 77], [74, 78], [75, 79], [76, 80]];



W = [V[:,k] for k=1:80]
hulls = [[[k for k=1+8h:8(h+1)]] for h=0:9]
cubes = AA(MKPOL)(DISTL([ W, hulls ]))

VIEW(STRUCT(cubes))

tangle = STRUCT(EXPLODER(cubes, 5,5,5))
VIEW(tangle)
