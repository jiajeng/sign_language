# sign_langu_TW
## content
- [Foder define(in nas)](#FOLDER) 
- [Experiment](#EXPERIMENT)
  - stimuli
  - process
  - event
 - [Data process](#DataProcess)

### FOLDER
- Experiment
  - E-prime (實驗程序-使用Eprime)
  - Event (實驗程序-甚麼時候開始某張投影片或是音檔)
- DATA
  - rawdata (.IMA file)
    - Interpreter (手語翻譯員, 專業人員)
      - SUBxxx
        - AO (A story oral)
        - AT (A story sign Language)
        - BO (B story oral)
        - BT (B story sign Language)
        - LOCALIZER (sentence)
        - REST (rest)
        - T1 (structual data)
    - Learner (手語學習者, 非專業人員)
      - SUBxxx
        - AO
        - AT
        - BO
        - BT
        - LOCALIZER
        - REST
        - T1
- code

### EXPERIMENT
#### STIMULI
- story(文章)  
  - B北投之旅(8段), 8分05秒
    fmri : 254 scan  
    *.\Experiment\E-prime\TSL Story 刺激材料\北投之旅_1113.docx*
  - A長壽的秘密(7段), 8分
    fmri : 242 scan  
    *.\Experiment\E-prime\TSL Story 刺激材料\長壽的秘密_1113.docx*
- localizer(語句)  
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*.\Experiment\E-prime\TSL localizer 刺激材料\localizer_CH_for_biling_project_2017oct.pdf*
![image](https://github.com/user-attachments/assets/c0f34188-0d41-44ed-9e63-4616e07a1de0)
- 影片檔放在E-prime程序資料夾(e.x .\Experiment\E-prime\TSL Story 刺激材料\灰色背景\Story fMRI E-prime\stim)裡的stim資料夾內。
#### PROCESS
![image](https://github.com/user-attachments/assets/870ef6b3-aae9-4b12-a1e3-837940417f35)
![image](https://github.com/user-attachments/assets/f0daa45f-fc09-4439-8154-f9fc8b768141)
&nbsp;&nbsp;&nbsp;&nbsp;Localizer  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;把所有句子都放在一起，包含口語跟手語的，總共48個影片檔，排序位置在EVENT檔案裡。

#### EVENT
沒有需要做反應，所以每個受試者的EVENT都一樣  
&nbsp;&nbsp;&nbsp;&nbsp; *.\Experiment\Event\fMRI_Localizer_scan數值.xlsx*     
&nbsp;&nbsp;&nbsp;&nbsp; *.\Experiment\Event\fMRI_Story_Scan參數.xlsx*   
放Eprime的程序，輸出的是甚麼，開始以及結束時間，以及對應的scan是多少。  
story 應該就不需要有event了 拿整段的資料就行  
localizer 才需要event  

### DataProcess

- step1  convert .IMA to .nii file (save in nas, no save in local)
   `save in LabDatajengsign_langa_TWDataniidata`  
    
- step2  get needed .nii file to local  
      
- step3  preprocessing image using conn
  `save in ./LabData/jeng/sign_langa_TW/conn_prj`
    
  每個condition分開跑一個project，project資料夾內有conn_xx_info.xlsx，表示跑了那些資料以及跑什麼流程     
     
- step4  calculate inter-subject correlation for each voxel  
  `save in ./LabData/jeng/sign_langa_TW/Data/ISCdataLearner(or Interpreter)`
    
  拿preprocessing完的資料(dswau.nii)，去跑每兩兩受試者的每個voxel的相關性，有四個condition(AO,AT,BO,BT)。     
  `folders...`    
  `.pairs_Z -- every subject pair fisher z(calculate from r, tanh(r)) `  
  `.pairs_corr -- every subject pair r value`  
  `.pairs -- empty`  
  `.group -- pairs mean and median data`
    
- step5  statistical
  `save in ./LabData/jeng/sign_langa_TW/Data/ISCdata...`  
    
  看condition之間有沒有差異  
  - 拿某個condition的某個voxel的資料來說，會把correlation matrix的某些pairs的r值乘上-1 (打亂數值)  
  - 把每個condition的每個voxel都隨機乘上-1，再對每個矩陣去做hyothesis的相減(or 相加) `假設有4個condition被納入，AO, BO, AT, BT，hypothesis是1, 1, -1, -1，oral-sign`
  - 取r值相減後的fisher z
  - 重複2000次，得到隨機的一個分布，再來看真實的數值位於這個分布的哪裡來決定說，這個hypothesis有沒有被接受
       
  計算condition之間的差異，有幾個hypothesis    
  - pro-learn for sign -- lang. pro(AT+BT)-learn(AT+BT) `.stats_pro-learn_sign`    
  - pro-learn for oral -- pro(AO+BO)-learn(AO+BO) `.stats_pro-learn_oral`    
  - oral-sign for pro -- pro  (AO+BO)-(AT+BT) `.stats_pro_oral-sign`    
  - oral-sign for learn -- learn  (AO+BO)-(AT+BT) `.stats_learn_oral-sign`
       
 > [!Note]
 > 如果不確定hypothesis是甚麼(哪個條件減哪個條件)，可以再資料夾內的threshypothesis.txt看   
 > 結果有用mean(ISC_hypothesis01_mean_p.nii)跟median(ISC_hypothesis01_median_p.nii)跑的，看mean就可以了。   
 > p是left tail(A-B), lmp是right tail(B-A)。   
   
  
