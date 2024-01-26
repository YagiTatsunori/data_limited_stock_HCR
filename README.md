Data-limited資源用の新しいHCRの開発用レポジトリ

## MSE_compare.R
reference_points.R, scenario_before_management.R, management_HCR.Rもダウンロード  
scenario_before_management_func("ICES" or "Japan",  
                                "one-way" or "roller-coaster" or "random",  
                                0.75 or 0.5 or 0.25,  
                                0.75 or 0.5 or 0.25)
1. 管理前のシナリオは、ICESか日本か設定
2. 1でICESを選んだ場合、どのシナリオを選ぶか設定
3. 1で日本を選んだ場合、初年度の資源量がB0の何倍かを設定
4. 1で日本を選んだ場合、最終年の資源量がB0の何倍かを設定

## RP_scenario_ICES.R
1. Fischerの論文に載っていたパラメータをもとに100年間、さまざまな漁獲圧で漁獲  
  
2. そのデータから、reference pointsを求める  
  
3. ICESの管理前シナリオ(one-way、roller-coaster)で個体群動態シミュレーション  
  
4. ICESのHCRで100年間の漁獲管理
