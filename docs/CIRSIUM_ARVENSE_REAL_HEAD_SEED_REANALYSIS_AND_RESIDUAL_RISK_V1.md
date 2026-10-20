# Q2 実データ再解析：送粉者排除で「推定種子食害」と「実際に損傷と記録した種子」は逆方向

**2026-10-20 | aza3 PR #50 | source-authenticated real *Cirsium arvense* 2022 field experiment reanalysis.**

## 元データを入手した：前の「Zenodo未入手」を更新

Karpate, Blažek, Kollross, Rech, Sam & Jorge (2025), *Pollinators, not predators, affect florivore survival in Cirsium arvense (Asteraceae) reveals an experiment based on Structural Causal Modelling*, Research Square **preprint**, DOI [10.21203/rs.3.rs-6634172/v1](https://doi.org/10.21203/rs.3.rs-6634172/v1)。公開著者アーカイブ [Zenodo 10.5281/zenodo.15364152](https://doi.org/10.5281/zenodo.15364152) は、実際には公開 **51,505-byte ZIP** （MD5 `eeab80e4c988891e10fb1cffd4c9c057`）、[著者GitHub `Yogita25/Multitrophic_Interactions`, `First_release`](https://github.com/Yogita25/Multitrophic_Interactions/tree/8a151c62ed069f775a16b9bc5fbdce209cbb8691) を関連ソースとして登録している。CC BY 4.0。アーカイブの発行日は2025年5月8日。前段の `probe_cirsium_arvense_2025_zenodo_v1.py` は当時metadataを回収できず、その後に今回 TinyFish で直接公式APIを解決し、**著者の公開GitHubから生の field CSV/R scriptsを取得できた**。

| オリジナルファイル | Git blob SHA-1（真正元ファイル認証） | 実際の行数 |
|---|---|---:|
| `seeds.csv` | `ca99c762162d15f006f0c73b4603ca71fd320a29` | **615頭花** |
| `Herbivore.csv` | `b0578e8c814b371e910da119c66e45575bc8b675` | **181昆虫観察行** |
| `Pollinators.csv` | `1723e96b8fc8ea755832daa786c387c5f3c188a8` | **120観察区間** |
| `Herbivores.R` | `6016c74ef3a8783a376f6a6f5e50a4c35f81a75a` | 著者による集計・モデル |

**元原稿は609成熟頭花と記すが、公開seeds.csvは615行**。著者のさらなるフィルタ・重複基準・ID不明など可能性があるため、なぜ6差が生じるかの原因は未確認。615が正しい生物学的最終Nとも、609が誤記とも決めない。地点・処理のある有効な行は612であり、`Site_Code × Block_Code × Treatment_Code`から推定できる植物単位164、サイト×ブロック45、対照とPollinator Exclusionが両方観測されたペア39。これは著者本文にある196個体と矛盾の可能性があるが、**元CSV由来の推定キーの対象範囲と本文の登録個体の範囲が違い得る**ため、同一の母集団と見なさない。頭花ごとの `Vial_Code` は反復し、現行の公開CSVだけでは独立頭花の物理標本IDへ一対一に接続できない。

## 非自明な会計上の恒等式：著者の「Eaten」は目視損傷カウントと独立ではない

全615行で、

```
Flowers = Viable_Seeds + NViable_Seeds + Damaged_Seeds + Missing_Seeds
Eaten_Seeds = Damaged_Seeds + Missing_Seeds          # actual author R code
            = Flowers - Viable_Seeds - NViable_Seeds # exact algebra
```

を確認。従って `Eaten_Seeds` は**観察された食痕だけではなく、種子・小花の会計残差を含み**、`Damaged`はこの合計から相殺される。**その数字だけから、受粉済み種子を幼虫が実際に食べた個数が増えたかを独立に判定できない**。

- 元ファイルで `Missing_Seeds<0` が **80/615**、`Eaten_Seeds<0` が **62/615**、`Eaten_Seeds=0` が **27/615**、`Eaten_Seeds>0` が **526/615**。
- オリジナル `Herbivores.R` は実際に `seeds.eaten <- seeds[seeds$Eaten_Seeds>0,]` としてから `Eaten`・昆虫数のモデルに入る。これは**正の推定食害だけを条件づけた解析**であって、事前登録全頭花の食害率を自動的に推定するものではない。
- `Missing<0` を0に丸めて偽の生物値を作らない。元の符号付き `Eaten` は **診断的計算**として残すだけで、確定した食害確率には使わない。

## 植物を独立単位としたサイト・ブロック内の実データ対照

4処理を4サイトで実施した2022年の野外実験から、種子の数えられた612頭花を、元の `Site,Block,Treatment` ごとに植物単位へ集約。**39のサイト×ブロックでのみ**、同じブロックにあるPollinator Exclusion（organza袋）と無処理Controlをペアで比較。各植物で `Σ seed outcome / Σ Flowers` とした。差は Pollinator Exclusion − Control（単位：割合の絶対値）。

| 各ブロック内の平均差 | 実データでの方向 | 断定できる意味 |
|---|---:|---|
| **`Viable_Seeds / Flowers`** | **−0.5091** | 袋掛け処理群では有効種子として計数された割合が大きく減少 |
| **`NViable_Seeds / Flowers`** | **+0.3739** | 非生存カテゴリの割合が増加。ただし残余の負値がある頭花では一部カテゴリの計数品質問題を別途調査 |
| **`Damaged_Seeds / Flowers`** | **−0.0439** | 確認された目視損傷カテゴリは**減る方向** |
| **`(Damaged+Missing)/Flowers`** | **+0.1352** | 著者の**残差含む推定食害**は逆に**増える方向**。これを真の物理食害と断定不可 |
| **`Missing_Seeds / Flowers`** | **+0.1791** | 残差増加が上記の食害指標の増加を完全に説明する会計上の関係 |
| **`Insects/head`** | **−1.67** | 種子ファイルの羽化昆虫数カテゴリは減少。昆虫目レベル・頭花ごとの真の幼虫摂食量ではない |

**再解析単位・誤差の区別：** `seeds.csv` の実データはほぼ決定的な方向差を持つ一方、元著者の群間割り付け・サイト依存・ブロックの検出対象数を保って**各サイト内で39ブロックを復元抽出**する20,000回の層化bootstrapで記述的不確実性を別に計算している（対応する数値は原本に基づくCI artifact/結果JSONに記録）。source code SHAを検証し、元ファイルの不完全分類・負の残差を公開のまま明示。これは**花粉媒介を独立に操作して推定した潜在アウトカム**ではない。

## 独立実装で計算したブロック層化bootstrapと、相手種の結合不能数

著者の4原本を直接取得した独立JS計算では、`39`の地点×ブロック対応ペアを各地点ごとに復元抽出して**20,000回**のノンパラメトリックbootstrapを実行した。括弧は95%パーセンタイル区間（**統計学的に登録された因果効果の信頼区間ではない**）。

| 相手排除袋掛け－通常開放の差 | 平均 | サイト内ブロック復元抽出95%区間 | 正／負のブロック |
|---|---:|---:|---|
| `Viable_Seeds / Flowers` | **−0.5091** | **[−0.5720, −0.4456]** | 0正・39負 |
| `NViable_Seeds / Flowers` | +0.3739 | [+0.2784, +0.4732] | 33正・6負 |
| **`Damaged_Seeds / Flowers`（目視カテゴリ）** | **−0.0439** | **[−0.0799, −0.0115]** | 13正・21負・5同 |
| **`(Damaged + Missing) / Flowers`（符号付き計数残余）** | **+0.1352** | **[+0.0287, +0.2311]** | 26正・13負 |
| `Missing_Seeds / Flowers` | +0.1791 | [+0.0841, +0.2688] | 31正・8負 |
| `Insects/head` | −1.6678 | [−2.7007, −0.8076] | 3正・28負・8同 |
| `Eaten`正の頭花に限った食害率 | +0.1052 | [+0.0222, +0.1883] | 28正・11負 |
| **`Missing >= 0` の計数整合部分集合（37対応ブロック）** | +0.1409 | [+0.0648, +0.2178] | 27正・10負 |

この区間は**最初の独立JS集計**から計算した。GitHub Actions上で実行する独立Python実装は元Git blob SHAを検証し、同じサイト・ブロック集計を再現するようにした。元の混合効果モデルや検定p値を再現したわけではない。負の `Missing` を除外する感度解析も非ランダムな事後選別で、真の食害率を特定するものではない。

**地理的に一地点の偶然だけではない：** 39対応ブロックを4地点へ分けても、平均の `Damaged_Seeds/Flowers` は全地点で負、残余 `Eaten_Seeds/Flowers` は全地点で正だった。対応ペア数は campus12、construction10、ditch9、rapseed8。それぞれの差は、`Damaged` が −0.0382／−0.0250／−0.1033／−0.0095、残余 `Eaten` が +0.0935／+0.2220／+0.1278／+0.0976。**4地点で数値の方向が一貫することは、本物の虫食害が逆転した証拠ではない。** 共通の種子会計定義による同じ制約も4地点へ繰り返し持ち込まれている。

**最重要：** 元の昆虫観察181行のうち `site × block × treatment × vial` で、種子台帳の単一頭花行に一意につながるのは**99行**だけ、**75行は同じvialの複数頭花が候補**、**7行は元台帳の当該キーに一致しない**。植物×相手ギルド×頭花の直接的な相互作用を検証できると思って181行を無条件に頭花へ結ぶと偽の共起を生む。さらに615元頭花から、有効なsite/treatment行だけでも**同一site×block×treatment×vialの重複キー134種**がある。

これがQ2で不足している要素を定量化する：観察・操作の研究は実在し、生物相手による実際の効果の方向も**原著の操作区別**としてある。しかし著者公開データだけでは相手種の正確な頭花内アクセスと最終成熟種子を**同じ生物学的個体IDで完全に接続できない**。従って`BITA+SCH+BALANCE+SLK` の「同一形態`z`について、相手の違いだけによる正味利得の符号反転」とは呼ばない。

## 生物学的な大事な帰結：送り手が減って見かけの被害率が上がることは防御の失敗ではない

**花粉媒介者を排除する網袋は、実際の送粉だけでなく、種子食者の産卵管による接近、寄生蜂の入り口、花の局所温湿度、観察可能性を同時に変え得る。** 原著で「送粉者がいると種子食者の発育成功が高まる」という生態仮説は合理的だが、この原データだけでは**花粉を通じた種子資源増加**と**網袋の物理的侵入制限**を区別できない。

今回の本当の実データから、相互作用研究で検証すべき機構は、

`pollinator access → fertilized achene resources`
と
`physical mesh/phyllary/spine entrance → pre-egg access`
を**別々に操作・測定**した上で、`true seed-feeder oviposition/larval feeding`と`plant's absolute filled and germinating achenes`を接続することである。

同じ物理的頭花構造 `z` が、送粉者にとって有利でも種子食者を誘引したり、その資源供給を助けたりするなら、構造を単純な「防御」や「誘引」と二分できない。けれどこれを Nature 級の先行未知の新法則と呼ぶためには、**実際の植物の構造z × 独立に識別した相手Pによる繁殖成功の反応面**を測る必要があり、C. arvense の袋処理の結果をそこへ格上げしない。

## データ統合上のさらに重要な制約

`Herbivore.csv` には実際の昆虫家族・個体数・観察段階があるが、種子CSVの同じ `Vial_Code` は同一植物内で複数の頭花に繰り返し出現する。**（地点×ブロック×処理×Vial）により、同じ昆虫と頭花の個別対応が一意とは限らない**。無条件にCSV結合して181昆虫を頭花615個へ1対1割り当てると昆虫数を二重計上し、頭花形態×相手ギルド効果の偽の証拠になる。解析コードは各昆虫行のリンクを一意・多候補・一致なしに分ける。

`Pollinators.csv` の120観察区間は主に対照・捕食者操作群で、送粉者排除の録画自体が同じ方法で収集されたわけではない。ゼロ訪花を観察したと偽装できない。真の頭花向き／刺長／総苞隙間はこのデータにない。

## GitHub 反映と実証上限

- [本物の元ファイルを毎回取得する再解析コード](../analysis/reanalyse_cirsium_arvense_2022_real_head_seed_bookkeeping_v1.py)
- `.github/workflows/reanalyse-cirsium-arvense-2022-real-seeds-v1.yml` は、著者の公開GitHubの**固定SHA 8a151c62...**から `seeds.csv`, `Herbivore.csv`, `Pollinators.csv`, `Herbivores.R` を取り直し、Git blob SHA-1全4件で真正性を判定して数値の再現性を確認する。
- [真正データからの数値・bootstrapと負の値の出所を記録したJSON](../data/evidence/cirsium_arvense_2022_seed_residual_and_partner_exclusion_result_v1.json) — CI出力を確認してから固定する。
- 2025年原著の査読誌掲載は未確認であり、**これは2025年の公刊前プレプリントの原データを2026年に独立再解析した結果**と表現する。

**実測として支持される結果：** 送粉者排除を伴う袋掛けでは種子の成立数も羽化昆虫も大きく減少し、計数上の残余 `Eaten` と目視 `Damaged` が逆に変化した。**未識別：** 昆虫の相手種が同じ頭花構造の利益を符号反転させたか、送粉が実際に食害個体の餌資源を増やしたか、頭花のレゴモジュールが適応進化したか、表現型消失後にDNA機能資源が残るか。三つの中心質問の②に有効な**一次解析の独立な負荷試験**として扱う。
