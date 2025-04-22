## 0830
| | spoken | sign |
|--|--|--|
|Interperter| | |
|Learner| | |

先分開story，Interpreter spoken跟sign的差異以及Learner spoken跟sign的差異  橫向 -->  
在來比spoken中Interpreter跟Learner的差異 以及sign中Interpreter跟Learner的差異 直向  

## 0902
發現Learner中的subject 18、19是一樣的，所以在step4ISC的程式中，在算corrcoef的地方加一個if，如果R值有等於1的話就是錯誤的，輸出是哪兩個一樣。

## 0903 
統計中，
原先只有算right tail(標為_p_), left tail(標為_lmp_), 做完permutation後只有跑right tail 的fdr沒有跑left tail的fdr，  
改成算兩邊的fdr-->/thres/的資料夾內會有mean_p跟mean_lmp的.nii檔。  
使用xjview來畫結果圖，算完fdr的時候順便把圖給話出來，並且存在result裡。

## 0904
出來的結果中都是left tail有值，也就是說condi1 z值是負數 condi2 z值是正數 的時候才有顯著差異性，
- story A , oral , Interpreter-Learner
![image](https://github.com/user-attachments/assets/41199acc-1f27-4976-a3d9-0b422344c6a2)
Interpreter的z是負數的，Learner的z是正數的，
