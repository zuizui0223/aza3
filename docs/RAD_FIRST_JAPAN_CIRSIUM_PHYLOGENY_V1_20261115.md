# RAD-first 日本産アザミ属・核DNA系統推定計画 v1

**2026-11-15 | 方法変更の提案 | 実試料・採集・発注は未認可**

## 新しい主線

ユーザーの優先方針を反映し、**新規取得する核DNA系統データの主役をRAD-seq（ddRADなどの制限酵素ベースの共通プロトコル）に置く**。

`全日本のRAD配列 → 共有相同座位を評価した核DNA樹／ネットワーク → 個体・集団構造 → 頭花形質の祖先状態／反復変化／連動`

既存Moreyra et al. (2025) 350核座位／Compositae1061系統は **external benchmark**。新規target captureを第一選択にしない。ただしRADの共有座位・倍数性QCで広域樹が成立しなければ、誤った一本樹を強制せず、RADの推論範囲を近縁クレードに限定し、既存核系統で広域を記述する。この切り替え条件をpilot前に固定する。

**注意:** 現行GROWの向きO×総苞構造Bの博士必達範囲を拡大したという意味ではない。遺伝データ取得法の優先順位をRAD-firstに変更する提案であり、GROW必達と全国128種の広域拡張は別の予算・時程ゲートで扱う。

## 1. 3層の採集パネル（計画値、確保済み試料数ではない）

| 層 | 新RADの目的 | 最初の案 | 開放条件 |
|---|---|---|---|
| Pilot | 異なる系統・2x/4x/6xで共通座位が得られるか | **24–48試料**。倍数性／地理／保存法をまたぐ代表試料＋技術反復 | 試料受領、許可、DNA QC、cytotype確認 |
| Species breadth | 日本の放散内で種間の主要分岐と頭花形質変化を推定 | 最初は既存38 paper conceptsに対応する**2–3個体/分類群程度**の感度パネルを候補化。推定に成功後に未採集分類群へ拡張 | Pilotで共有座位、倍数性・系統樹頑健性を確認 |
| Focal population depth | 個体多型・introgression・祖先多型・形質の組替え | 形質差が検証された2–4対象系の2–3集団×約8–15個体/集団を**候補**とする | 現地実在、権利、同一個体の形態記録と採集可能性 |

全国128種・228代表資料の候補は[未マージPR #3](https://github.com/zuizui0223/aza3/pull/3)の拡張構想であり、RAD実測も試料確保も意味しない。検証前に`128 × 3`などを確定発注しない。また先行 `core190` / `298` は現在のpilot採集数ではない。

最初に優先する対照は、異なる倍数性で既知染色体情報がある分類群、向き・総苞の変化枝に近い集団、同種内で複数形質が共存する集団。野外試料の実個体倍数性は文献上の種名から一律に代入しない。

## 2. 具体的なRADライブラリ工程

1. 対象個体ごとにDNA用の葉を確保し、独立にvoucher／写真／O・B・T・S形質を記録。
2. 適切な生葉があれば同じ個体／各集団でフローサイトメトリーを実施し、必要な基準染色体数と組み合わせてcytotypeを分類。**2C DNA量≠倍数性の自動確定**。
3. DNA品質（量、分解、汚染）を判定後、候補制限酵素の二重消化 → barcode／adapter ligation → サイズ選別 → PCR／精製 → pool → QC。酵素・fragment window・目標深度はpilotで選ぶ。ライブラリとシーケンサーの仕様は担当実験施設と決める。
4. 試料の系統・倍数性・形質が調製日・index・laneと完全に重ならないよう割付。独立DNA抽出／技術反復をpilotに含める。
5. 原FASTQ、sample sheet、酵素、サイズ選別条件、調製ロット、読長、QC、データ利用条件を記録。

RADは必ずしも参照ゲノムを要求しない。ゲノムサイズ・倍数性により最適酵素条件と座位数が変わる可能性をpilotで検証する。ライブラリを別の制限酵素系で作ると共通座位を直接統合しにくいので、成功後は同じ設計を維持する。

## 3. 解析の順序（「RADを主樹にする」が成立する条件）

### A. 同一RADセットを処理

FASTQ → demultiplex/QC → de novo loci clustering（例：ipyrad、パラメータ感度必須）→ 相同座位候補 → SNP／locus配列セット。

- locus occupancy／個体の欠損率／深度分布を、種・クレード・cytotype別に出す。
- 制限サイトの変異によるdropout、極端な深度、過剰ヘテロ接合、paralog/homeolog collapseを検査。
- 欠損閾値を1通りで固定せず、**緩い／中程度／厳しい共有座位セット**でトポロジー感度を検査（ただし結果を見て閾値を都合よく変更しない）。

### B. 2xの核DNA系統を先に安定化

まず同祖先と判断できる単一コピーRAD lociを核として二倍体パネルの複数樹を推定する。

- full RAD loci alignment と SNP-only（ascertainment／linkage考慮）を混同しない。
- 連結ML樹（例：IQ-TREE 2）と quartet/共祖先法の感度を確認。短いRAD locusごとのgene treeは情報不足の可能性があるので、ASTRALを機械的に正解としない。
- 樹形に対する系統的に偏った欠損、集団サンプリング、外群選定の影響を点検。

### C. 4x/6xの参加条件を別にする

- フローサイトメトリー・染色体情報と、二倍体／多倍体の祖先関係を照合。
- 個体の倍数性に応じてアレル数・用量と不確実性を取り込むgenotype calling（例：polyRAD）。**polyRADは倍数性を推定しないし、homeologの正確な分離も保証しない**。
- 四倍体AAABのような状態を二倍体ABと同じ確定genotypeとして統合しない。
- 異質倍数体で祖先由来コピーを識別できない時、単一terminalによる「確定種樹」を強制しない。二倍体骨格＋polyploid placement／ネットワーク／サブゲノム別歴史を報告。
- 2x-only / 同一cytotype-only / mixed strict ortholog で主要分岐が変化しないか比較。

### D. RAD由来の系統樹の主張上限を判定

| 結果 | 許される結論 |
|---|---|
| 広域で十分な共有相同座位、樹形が欠損・倍数性・手法に頑健 | **RAD-first の日本産核DNA主樹** |
| 近縁クレードでは頑健、遠縁間では共有座位が乏しい | **クレード別RAD主樹＋広域は未同定** |
| polyploidが混合祖先を示す／homeologが分離できない | **2x主樹＋polyploid網状系統／保留** |
| RAD樹と既存350核遺伝子樹が不一致 | **座位・欠損・ILS・交雑・倍数性を比較**。RADを無条件優先しない |

**枝長と年代:** RADのsubstitutions/siteやSNP-only由来の長さをそのままMyrと読まない。年代推定と形質進化速度は、年代較正・座位sampling・ascertainmentを検証してから行う。

## 4. 葉緑体DNAは別の母系史

RADの核DNA結果と、同一／対応個体の葉緑体haplotypeを**別々に**推定し、矛盾を記録する。

- 同質倍数化のみで母系系統と核系統が自動的に矛盾するわけではない。
- 異質倍数体の母系、交雑後のchloroplast capture、ILS、葉緑体の低分解能が主な候補。
- 核対葉緑体の不一致だけではgene flow/倍数化を確定せず、RADの遺伝子流動／集団構造とcytotypeデータを照合。
- 葉緑体をRADの核SNPに連結し「一つのDNA系統樹」と称しない。

## 5. 形質の問いへ戻す

同一個体の実測を最優先：O（重力基準角、蕾・開花・結実時期）、B（総苞開度・間隙・被覆）、T（真正トゲ長・角度・硬さ）、S（粘着力／腺／分泌物）、花色・撮影情報を保存。

RAD核系統の樹形不確実性と多型/欠測を伝播して：
- 祖先状態・最低変化数（独立起源数とは区別）；
- 相互の変化枝の重複と非同期性；
- 分岐史・集団構造・遺伝子流動でどこまで形質差が説明されるか；
- 祖先多型／交雑／別々の由来の競合仮説。

RADの形質連関だけでは選択、原因遺伝子、遺伝的モジュールの独立性、機能・適応度を主張しない。機構とfitnessは別に実験する。

## 今の実務ゲート

[2026-11-15の物理試料受領監査](RADSEQ_MATERIAL_READINESS_GATE_20261115_JA.md)は公開台帳3件すべて **実個体登録0件**。研究室や協力者に物理試料が存在しないという意味ではない。

1. まず所有／受領／使用許可／voucher／葉の保存状態／個体と形質の対応を非公開台帳で確定。
2. 得られた試料から **24–48候補**を系統・倍数性をまたいで選び、RAD切断・ライブラリとde novo clusteringの予備試験を実施。
3. Pilotの共通座位と品質が成立してから、全国の拡張パネルを定量的に設計。
4. 採集・利用・配列の許可は本書の作成や承認で自動的に与えられない。

## 方法上の参考文献

- Karbstein et al. (2024) [A roadmap of phylogenomic methods for studying polyploid plant genera](https://pmc.ncbi.nlm.nih.gov/articles/PMC11342234/) — RADの種間限界、異質倍数体、homeolog。
- Clark et al.; [polyRAD R package](https://cran.rstudio.org/web/packages/polyRAD/index.html) — ploidy-aware allele-dosage genotype uncertainty。
- [ipyrad reference](https://ipyrad.readthedocs.io/en/master/6-params.html) — de novo mode without reference。
- [Data-driven SNP phylogenomics guidelines](https://pmc.ncbi.nlm.nih.gov/articles/PMC11610416/) — locus filtering、branch lengthと祖先形質。
- Moreyra et al. (2025) — East Asian Cirsium nuclear backbone as independent comparison.

**意思決定:** RAD-firstを研究データ取得の主線にする。RADで推定不能な広域・異質倍数体の祖先関係は、推定可能と装わず別階層に保留する。
