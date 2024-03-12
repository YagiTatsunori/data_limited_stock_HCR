Data-limited資源（２系資源）用の新しいHCRの開発用レポジトリ

## functions.R
stock_parameters()でインプットデータから、年齢ごとの体長などの必要なパラメータを計算
reference_points()で、管理基準値を計算
MSE_result()で、管理前の資源動態と、HCRによる管理を実行
performance_MSE()で、MSEの結果のパフォーマンス計算を実行

  
MSE_result(scenario_organization = "ICES" or "Japan", # 管理前のシナリオは、ICESか日本か設定
           scenario = "one-way" or "roller-coaster" or "random", # ICESを選んだ場合、どのシナリオを選ぶか設定
           start = 0.75 or 0.5 or 0.25, # 日本を選んだ場合、初年度の資源量がB0の何倍かを設定
           end = 0.75 or 0.5 or 0.25, # 日本を選んだ場合、最終年の資源量がB0の何倍かを設定
           sd_r = 0.6, # 再生産誤差
           sd_i = 0.2, # 資源量指標値の観測誤差
           sd_l = 0.1, # 毎年の体長の平均値の誤差
           sim = 1000) # シミュレーションの回数

## bio_param.R
von Bertalanffyの成長式でのkの値によってICESのHCRはパラメータを調節する
なのでさまざまなkの値を持つ、以下の資源のパラメータを設定する
pollack (Pollachius pollachius; k = 0.19)
Thornback ray (Raja clavata; k = 0.09)
Brill (Scophthalmus rhombus; k = 0.38)
Plaice (Pleuronectes platessa; k = 0.23)
