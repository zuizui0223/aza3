# 三つの進化的質問をアザミで識別する：記憶 × 利得 × 再到達可能性

**2026-10-19 | 研究プログラムの新規統合案 | 所管：aza3（Cirsium 実証） | 状態：prospective。未許可の採集・操作は実施しない。**

## 結論（3問を一つに統合するが、3問の証明を混同しない）

統一的な問いは**「進化史によって保存・交換された繁殖形質の構築資源が、現在の昆虫群集で機能的利益をもち、失われた表現型を再形成・維持できる条件は何か」**。

同じものとして扱ってはいけない3層：

1. **記憶（Q1）**：過去の進化史と形質が、現在の頭花形質の似かたをどの程度説明し、その程度は別系統に予測として移植可能か。
2. **利得（Q2）**：実際の構造・色・匂い・姿勢が、それぞれの相互作用相手と環境下で、有効受粉・種子食・寄生蜂のアクセス・充実痩果のどれを変えるか。
3. **再到達可能性（Q3）**：外見で形質が消えても、遺伝的・発生的資源が同一個体・集団に保持され、発現可能で、集団内で遺伝・持続するか。

`past dispersal / ancestry H -> functional resource inventory G -> regulatory/developmental compatibility F(G,E) -> real capitulum phenotype Z -> guild-specific access A(Z,P) -> reproductive outcome W -> descendants`

枝分かれが重要：生態環境・訪花群集 `P` は来訪頻度と `W` を変えるだけでなく、将来の遺伝子頻度を選択する可能性がある。遺伝子流動 `J` は再出現に見える新しい資源を輸入し得る。**P の変化を過去の G の喪失・復活と自動的に結びつけない**。

## 現在の最小実証ゲート：JPN36 *Cirsium sieboldii*（EAzamiの原本から検算）

[EAzami JPN36の現場・ゲノム橋渡し判定](THREE_QUESTIONS_JPN36_MINIMUM_EMPIRICAL_MILESTONE_V1.md) が、3問のうち最初の **Q2実在機構**への有力候補である。日本若年放散の末端総苞姿勢イベント（UFBoot局在0.754）は形態史の背景。実測ギルド機能はまだない。EAzamiの凍結済み計画は非操作基準 **30独立個体**、安全性確認各処理6頭花以上、許可後だけ12対/24頭花の操作試験。頭花操作が計数開花個体の2%以内なので、**24頭花全体を対象にする必要条件は少なくとも1,200開花個体**。ただし1,200は個体数上の条件のみで、サイト別保全・権利・安全性の承認は現在**未取得**。

aza3の *C. sieboldii* Phase-A **40予定個体**は別の分子/同一個体高次LEGO計画であり、EAzamiの24操作個体や30基準個体と実際に重複すると確証されたIDは**0**。両計画の数を合算して実試料とせず、同一植物×頭花×実験日×分子採取の明示ID照合が得られるまで Q2×Q3 の統合実証には昇格しない。`analysis/audit_jpn36_sieboldii_three_question_bridge_feasibility_v1.py` は公開固定原本に対する源データ＋保全ゲートの再現検証を実施する。実際の生物学的効果量を推定するものではない。

**何が実証可能かの優先順位**：基準観察・安全性・許可がそろえば **Q2 総苞の操作可能性**を先に評価。**Q1 別系統への記憶移植**はこの1種では不可能、**Q3 表現型消失後の再進化**は遺伝的資源/祖先状態/再発現・遺伝を別途実証する必要がある。今すぐ3問が一つの野外操作で完了するという解釈を採用しない。

### 横断リポジトリの根拠と上限（それぞれの凍結済み結論を変更しない）

| 上流 | 支持されている情報 | aza3で検証すべき独立な帰結 |
|---|---|---|
| [turnover](https://github.com/zuizui0223/turnover) | 201形質×科システム、45植物科、12連続形質。形質記憶の科ごとの反復性 0.150（95%区間 0.020–0.286）、leave-one-family-out 形質名の予測改善 **−0.024**（S3）と **−0.032**（prune）。**同じ形質の記憶強度は別科へ単純移植できない**。一般化済み時間×空間4セルの主設計は source gate HOLD。 | アザミ内で「向き、総苞、粘着、花色」の記憶勾配を下位系統別に推定、**未使用の下位系統を丸ごと保留した予測**。属内結果で「科間移植」を証明したと呼ばない。 |
| [chun](https://github.com/zuizui0223/chun) | 花色の類似性は系統分岐とともにしばしば低下するが、最も情報を残す表現型の粒度は普遍的でない。可視色の類似性≠同じ分子経路。 | 頭花の連続色・色素化学・遺伝型をそれぞれ独立表現して、どの解像度で歴史を予測できるか。CamelliaをCirsiumとプールした独立データと見なさない。 |
| [EAzami](https://github.com/zuizui0223/EAzami) | 日本の若いCirsium放散で向き ML 最小6、総苞姿勢3、粘着5の状態変化、相対的深さが異なる。**同期枝位置3組すべて堅牢な共変化証拠なし**。現生気候相関は特定の起源時環境を復元できない。 | 反復した形態は候補だが、**遺伝資源の再利用や適応・可逆性はまだ未確認**。系統独立性・祖先状態不確実性を明示。 |
| [bita](https://github.com/zuizui0223/bita) | 実際の相互作用データでは、アクセスが制約されると正当利用だけ低下し、迂回・盗蜜は比較的残り得る。証拠は鳥花系のネットワークなどであり、Cirsiumへの効果量の移植は不可。 | 送粉者・種子食者・寄生蜂が**頭花の同じ物理インターフェース**を使うかと、失敗後に迂回するかを同定する。 |
| [sch](https://github.com/zuizui0223/sch) | **多機能性≠機能衝突**。同じ軸に対する相互作用相手別の選択勾配が本当に逆向きか、相手操作を交差した `W00(z)...W11(z)` が必要。 | 一つの実測 `z`（例：頭花の重力基準角度）で、花粉機能と種子食・降雨機能を分け、対立も**一致する最適値**も検証する。 |
| [balance](https://github.com/zuizui0223/balance) | 妥協領域の存続/履歴依存を区別。複数機能が同じ器官にあることは、共有最適値・最適値対立・独立軸を識別しない。 | 複数機能の実際の適応度座標・頭花発育期・同じ器官部品を特定するまでは BALANCE 状態と認定しない。 |
| [slk](https://github.com/zuizui0223/slk) | `L`（衝突）→`R`（分化で回収できる利得）→`Phi=R-K`（構造費用差）→局所到達可能性→侵入・固定は**別のゲート**。 | `L>0` だけで頭花の独立モジュールが有利とは呼ばない。形態・器官発生費用と遺伝的組立可能性を測る。 |
| [anto](https://github.com/zuizui0223/anto) | 花色が同じ白でも共同保持資源 `G` は異なり得る。**3資源の周辺・ペア頻度が同じでも、三つ共同保持率は0〜25%**。完全な現状態を測れば、理論モデルでは過去の付加予測力が消える；不完全測定では残る。すべて理論であり再発現の野外実証ではない。 | 可視表現型でなく**同じ個体内の機能資源の共同保持と発現・相互適合性**を評価。歴史の独立予測力を、分子状態の測定精度とともに評価。 |
| [egwe](https://github.com/zuizui0223/egwe) | 実装源に固定した256種子×4設定（計1,024**計算軌跡**）の理論実験では、初回発現と3世代存続・12世代占有・資源の完全喪失が異なる。低保有率では一過性再出現の後に失われやすい。**植物の実測ではない**。 | 白花で一回色素が誘導されても、遺伝する／複数世代で維持されるとは呼ばない。 |
| [izu-core](https://github.com/zuizui0223/izu-core) | Model 3で、進化順序の先行は必須の因果順序でない。初期のアレル支持集合、分散と再導入の有無が到達可能性を制約。Q2 384軌跡の2主要効果は `both_inconclusive`。 | 再獲得の祖先系譜・新突然変異・移入・組換えを分ける。花粉媒介者の復帰＝色素対立遺伝子の自動回復とはしない。 |
| [aza3](https://github.com/zuizui0223/aza3) | 現行実証所管。2026若い琉球・台湾系統の公刊核遺伝系統；花色・色素・同一頭花モジュールの prospective design、40個体 `C. sieboldii` 多モジュール、既定298個体 Phase-A、既存25葉RNA登録。近縁DFR/ANR・ANS/FLS候補の7比較は機能的酵素として**未確証**。 | 本計画の実測・同一個体ジョイン・染色体/ゲノム/遺伝子型/色素・頭花機能/充実痩果は**すべてaza3に置く**。 |

## 2026-10-19 実データで分かったQ1の接続範囲

詳細：[日本の頭花状態と東アジア候補を比較した源データ監査](CIRSIUM_HEAD_STATE_TRANSFER_VS_EVOLUTIONARY_MEMORY_20261019.md)。

- 日本 Moreyra38名簿のうち NMNS の頭花向き二状態を**厳密に同じ定義で**比較できたのは19分類群（上13・下6）。別に上向き/上昇1、索引と詳細記述の矛盾1、unknown3は昇格させない。
- 東アジアの探索非重複19候補から、向きが事前に解決した11種（上6・下5）を固定。日本19で学んだ「上向き比率13/19」の外部Brierは**0.267**、50:50の無情報比較では**0.250**。これは種の**状態頻度**の予測陰性対照であり、**進化的記憶の非移植性の検定ではない**。
- 公刊核S1標本258行と東アジア11種を照合して、核DNA標本記録が正式な名称・voucher・地域と整合するのは**2種だけ**（*C. argyracanthum* 中国・*C. vlassovianum* モンゴル）。両方とも植物誌では上向き。*C. lidjiangense* に見える第三の一致は**トルコ標本と提出名不一致**で棄却。残り8種はこのソースに正確な一致なし。
- **独立したGBIF座標品質ゲート**を重ねると、出現点を確認できる向き候補は下向き *C. schantarense*（8点）と上向き *C. vlassovianum*（8点）だけ。**核系統標本と分布点の両方に安全な一致があるのは上向きの後者1種だけ**で、下向きは0種。*C. schantarense* の核DNAが未知という意味ではない。**Bae 2015（DOI 10.5352/JLS.2015.25.2.243）**は下向き *C. schantarense* と上向き *C. japonicum* が比較的近い位置にある **linked ITS/nrDNA** 系統を既に報告している。次はその配列・voucher・種名を回収し、独立の多数核座位で同一枝長の系統に接続可能か調べる。分布点と独立rDNAは系統起源回数や形態記憶の強さではない。
- したがってQ1の最初の遺伝・系統ボトルネックは**下向き候補を含む、真に出典照合できる核系統の複数独立分枝と同一尺度の枝長**。現段階で外部系統の記憶 `rho`、転移利得や遺伝資源再利用を推定してはならない。核標本の植物誌上の種状態はDNA個体に直接測った重力角度でもない。

この情報は既存の生物学的実験の許可、`EAzami` の形質状態・系統結論、Q2の機能実験・Q3の分子実験の凍結済み判断を変更しない。

## 分類群の役割分担（すべてCirsium、ただし独立性を捏造しない）

| 検証目的 | 第一候補 | なぜ必要か | 明確な限界 |
|---|---|---|---|
| Q1 同じ形質の記憶が別系統へ運べるか | 日本側 `Nipponocirsium`、琉球 `Arenicola`、台湾 `Sinocirsium` に加え、別大系統の北米Cirsiumを**可能なら独立ホールドアウト** | 同定した相同形質を各系統で共通測定、phylogenetic patristic距離×形質非類似度を比較 | 2026年公刊の三亜節は共通祖先をもつ。属内の成功だけで植物科間一般性を主張できない。属内十分な種・変異がない形質は HOLD。 |
| Q2 機能と相互作用相手 | EAzami向き・総苞・粘着の状態差をもつ複数同属集団；既存 JPN36 は観察/feasibilityのみ | 同じ頭花ステージ・時刻で、異なる昆虫ギルドが実際にどこを通るかを見る | Pi 3クラス昆虫目検出は送粉者/寄生蜂/種子食者の役割同定には不足；真のトゲではない輪郭突出を代用不可。 |
| Q3 表現型消失後の到達可能性 | `C. brevicaule`白 × `C. irumtiense`淡紅、独立の `C. japonicum var. takaoense`白/青紫、`C. sieboldii`同一個体LEGO | 白・着色、発現、遺伝、祖先/移入、共同部品の区別 | `brevicaule→irumtiense` の実際の白花化→淡紅色再獲得**順序は未確定**。別種であることは既知で、それだけでは祖先帰りといえない。`takaoense` 既存白/青紫3+3葉RNAは場所と表現型が交絡。 |

**重要な新知見とリスク：** Chang et al. (2026) `10.1186/s12870-026-08097-6` に基づく37アクセッションの遺伝系統により、`C. brevicaule` と `C. irumtiense` は宮古海峡で分かれた別系統であり、2C DNA量は報告上約1.93 pgと約0.69 pgと大きく異なる。この比較は因果的な「元は白で戻ったピンク」の証拠ではない。島をまたぐ移住史・細胞型・遺伝子流動・祖先状態不確実性を独立に扱う。

## 実際の頭花データが教える、最初の「識別できる比較」

最新の原著 Chang et al. (2026, DOI `10.1186/s12870-026-08097-6`) は、名前・形態・色を区別する際に重要な**自然な陰性対照**を提供している。ただし、これらは著者の形態比較の平均値であり、私たちの追試や交絡調整済み効果ではない。

| 比較 | 原著の実測された頭花関係 | 役割 | してはいけない推論 |
|---|---|---|---|
| **var. takaoense 白 vs 青紫**（同一分類群） | 内総苞片／外総苞片の長さ比 **2.71（白）対2.85（青紫）**。原著は他の主要小花構造の有意な相違を報告していない | **色だけが異なり頭花構造が比較的似ている自然対照候補**。Q2で色の誘引と外形アクセスを分離するため優先して現地内比較を設計する | 有意差がない≠形態同等性の証明。同一集団内の白/紫ペアではなく、既存6葉RNAは標高・地理と色型が交絡するため色だけの因果比較ではない |
| **C. brevicaule 白 vs C. irumtiense 青紫**（異なる琉球種） | 頭花幅 **1.79対1.93 cm**、内／外総苞片比 **2.24対1.44**、小花長 **2.94対2.59 cm**。ゲノム2C量も **1.93対0.69 pg** | Q3で祖先系譜と共同資源の由来を分ける、**系統・形態・ゲノムの束として異なる比較** | 組合せ差を色の効果と呼ばない。この2種が直接祖先→派生の白→紫再獲得イベントを示したわけではない |
| **日本の C. sieboldii 複数形態** | Phase-A `40` は新規採集予定値で、実際の同一個体の色・向き・総苞・粘着・機能ゲノムはまだ同時測定していない | Q1×Q3で、見かけ上の相関と別々の遺伝資源の共同保持を区別する第一候補 | 40個体＝40系統独立起源、あるいは組合せが存在する証拠と呼ばない |

ここから**もっとも効率がよい実施順**は、(a) var. takaoense の白／青紫を**同一集団・同じ発育期**で直接測れるかを現地適格性だけ先に確認し、(b) `C. sieboldii` で色と形態の多軸変異が本当に**同じ個体群**に存在するか検証し、(c) 琉球2種のゲノム/花器官の発現の違いを系統史と切り離すことである。地理で完全分離している場合はそれ自体が厳格な実験的限界となるので、補償的に高標高BP・低標高Wの追加候補を事前に探す。**既存葉RNAの6 W/BPサンプルだけで色×昆虫の因果実験をしたとは呼ばない。**

**原著の書き方に関する監査事項：** 2026年論文は要旨/結論で var. takaoense 色差をアントシアニン発現・送粉者選好と結びつける強い表現を用いるが、今回source-lockedされた25試料は**若い葉**であり、そこから小花の酵素活性・花器官RNA・同一個体の昆虫選好実験を再現することはできない。原著の表現の強さと**入手済みの検証可能な一次観察**を区別すること。現行の新規アザミ実験は、その不足部分を明示的に測る。

## 三つのワークパッケージ（同じ個体を最大限共有）

### WP1 — 形質の記憶と他系統への予測移植

- **測定粒度**：頭花の重力基準角度、花頭長径/短径、総苞片の姿勢/真の刺長/刺の方向、粘着分泌、floret可視+UV色・HPLC/LC-MS 色素を別表現。発育期、頭花順位、花の性、撮影/計測誤差を標準化。方向は写真の上下から推定しない。
- **応答量**：各系統内の `rho_clade,trait = corr(phylogenetic distance, trait dissimilarity)` を、系統依存・非独立ペア・標本数・年代の不確実性に配慮して算出する。rhoをOU `alpha`・半減期・進化速度に読み替えない。
- **主要予測**：他のCirsium系統で学習した「同じ形質名」の記憶勾配を、未使用系統へ適用したとき、系統固有・形質×系統モデルより予測改善が乏しいか。定量予測は `leave-one-lineage-out` と holdoutされた植物個体群で評価。
- **反証と停止**：系統内の相同形質に十分な変動や複数の明瞭な遷移がなければ memory 推定不可。新大系統の樹木・種数・形質を見てから追加/除外を決めたモデルで普遍性を称さない。属外の同一測定を独立に確認しなければ **Q1 genus-level portability only** と記述。

### WP1現有系統の数学的ゲート：琉球2種だけでは「記憶の傾き」が出ない

2026年の公刊系統研究における `Arenicola` は **C. brevicaule／C. irumtiense の2種**である。種間組合せは `2×1/2=1` ペアに過ぎない。したがって、この2種間の頭花色・向きの違いがどれほど明瞭でも、**系統距離と形質差の相関（進化的記憶勾配）自体は推定不可能**である。1本の距離と1本の差には順位相関も信頼区間も定義できない。各島の集団・個体を増やしても、それを独立な**種間進化起源**に置き換えることはできない。

台湾 `Sinocirsium` の5系統、日本・台湾 `Nipponocirsium` の数系統も形質が十分変化しなければ記憶推定が安定しない。**最初に**同じ器官の測定と十分な系統分岐・形質変異の情報量を診断し、その後に独立な日本の種群と北米大系統などを固定ホールドアウトとして選び直す。支持不足なら **Q1 HOLD / Q2・Q3 続行** が正しい判定であり、Q1を捻り出すために同じ島の個体を進化的独立反復と数えない。

### WP2 — 相互作用相手が構造の利得を変えるか（最初に実証すべき主柱）

- **生態学的単位**：植生・開花資源・時刻・蕾/開花/結実段階・真の頭花順位と性を揃えた、複数独立母植物上の**全事前登録頭花**。連続録画の撮影分母には昆虫未到来の頭花も保持する。
- **昆虫が実際にすること**：正当な送粉者の柱頭接触＋花粉沈着、種子食者の外側接触/失敗/迂回＋真の産卵、寄生蜂の幼虫探索/産卵と種子の不可逆損傷との先後。GloBI関係やDiptera/Hymenoptera/Lepidoptera自動分類は、その種が頭花内で本当にした行動の代替にはならない。
- **機構**：世界方向 vs 頭花軸方向の侵入、頭花を構成する真の刺・総苞の隙間を校正して測る。飛来の誘引/香り・着地・入口利用・産卵を区別する。実際の侵入失敗の全試行分母を保持。初回阻害だけで完全防御とせず、**迂回後の最終アクセスと卵/花粉**を分ける。
- **主要な因果実験**：観察でギルドと物理接触が確認でき、現場の許可/検出力ゲートを通った後に、まず**安全に再現可能な一軸 `z` の3段階以上の変化 × 適切なsham、相手の存在/不在（または独立に操作可能な相手機能）**を実施。目標 `W_{00}(z),W_{10}(z),W_{01}(z),W_{11}(z)` は同じ頭花の最終**充実痩果数**を主要尺度とする。全頭花の適応度の介入解析では「昆虫が来た」「卵が羽化した」頭花だけを抽出しない。操作単位が植物なら**植物が独立反復**であって頭花数ではない。
- **別の有力競合**：rain/UVによる花粉保護、既授粉の匂い減少、頭花大小と総苞順位、他Cirsiumの同時開花、花色化学、実際の刺に触れないギルド。トゲ切除・粘着除去はそれ自体で組織損傷/匂いを変えるため、shamと機構監査なしでは実験しない。
- **SCH/BITA/SLK判定**：同じ形質で相手別の最適値が逆向きなら「衝突候補」、同じ向きなら「機能整合」。相手ごとの妥協した状態の表面と、純粋機能最適値を混同しない。損失L、回収可能利得R、構造費用Kを直接同定しなければ `Phi=R-K>0` とも独立モジュールの自然選択とも呼ばない。

### WP3 — 見えなくなった頭花色は、遺伝的に戻れるか

- **現在の白花は機能的ゼロではない**。色素構造遺伝子の完全ORF、転写調節MYB/bHLH/WD40、DFR/ANS/ANR/FLSの実際のパラログと酵素機能、色素中間体、花弁ではなく**アザミの小花**の時期別発現、必要なら非花組織の保持機能を区別する。
- **識別する5歴史**：祖先着色の保持／祖先の花色消失の後の内生的再発現／新規突然変異による再構築／色素アレルの浸透交雑／着色系統の単なる移住。さらに別色素系と一時的環境誘導は対照。見た目だけで白→ピンク再進化とはしない。
- **有効な一次測定**：同じ個体の小花+総苞の色素と発現、全体の ancestry/cytotype、候補座位のハプロタイプ/コピー数/局所gene tree、実際のタンパク質・酵素機能と共通環境。候補相同配列が見つかっただけ、葉RNAで発現しただけ、保存されたmotifだけは**発現可能性の検証に不十分**。2026-10-19時点のaza3 DFR/ANR–ANS/FLS の7比較はいずれも同定区間が0を含み機能確定不能。
- **実際の到達可能性と存続**：白花系統/集団が同じ背景・条件で色素経路を再活性化できるか（許可された共通園・刺激処理・遺伝学的相補・交配、陰性/陽性対照）、**花の色素化学が本当に変わったか**、その変化が遺伝/分離するか、昆虫相互作用を変えて充実痩果に利得があるか、次世代でも維持されるかを**別々に報告**する。
- **最高リスクのテスト**：単一時点の白花だけでは区別できない `G` の三資源共同保持に、独立機能アッセイを割り当てる。過去の系統/地理情報による次状態予測改善が**可視色のみでは残り、測定した機能資源を加えると消えるか**を保留集団で比較。消えなければ測定不足、進化史の情報、またはモデル誤差を検討し、非マルコフ性を直ちに認定しない。
- **停止**：分子・祖先・繁殖世代のいずれかが欠ければ `LATENT_CAPABILITY_ONLY`、`PIGMENT_INDUCIBLE_ONLY`、`HERITABLE_REGAIN_UNPROVEN` などの明示的な上限。別種の既存色素状態は「再出現」ではない。

## 3問の統合を成立させる「共通データ」対「追加データ」

**共有の基盤（各個体・各頭花ID）**：植物集団の位置、核/葉緑体 ancestry、cytotype、頭花の性/順位/発育期、連続形態5軸、floretの可視+UVと色素、小花/総苞の発現と化学、昆虫の同時録画/同定と行動、袋掛け/開放/shamを許可状況に応じて区別、充実痩果の最終評価。

**Q1専用**：複数の独立な形質起源・他下位系統の樹木/種単位ホールドアウト・祖先状態・形質相同性。

**Q2専用**：安全に操作できる共有形態軸、昆虫相手/降雨条件の交差割付、同一頭花での真の侵入試行＋産卵/花粉成功、個体当たりの絶対成熟痩果。

**Q3専用**：花器官に限った機能的遺伝資源（酵素/制御/遺伝）の保持・共同存在・家系と再組立可能性、実際の形質の遺伝的伝達・持続。

上記以外を一度に採集して網羅的解析と称さない。**異なる群集での同じ構造の利得逆転**と**失われた形質が戻ること**は、反応変数も操作も違う。

## 実行順と論文構成

**Phase 0（すぐ）**：既存 source trees（EAzami、公刊2026年Cirsium系統）の属内ホールドアウト適格性、同じ形態の測定単位を照合。v5遺伝子相同性の未確証状態を保持。生態/分子/年代を装った合成値を作らない。

**Phase 1（最小の現地実証）**：頭花本物の向き・総苞形態・色素と、ギルド特定の**自然接近→接触→産卵/花粉**を同一集団で確認。独立植物を検出力単位とし、有効動画努力分母と最終痩果数を揃える。観測上の正の経路がなければ複雑な実験へ進まない。

**Phase 2（主要生態学論文）**：Phase 1 の実在ルートに基づき、許可済み単一主形態軸×相互作用相手の計画を登録し、費用/利益が**同じ軸で本当に対立・一致・迂回するか**を試験。各ギルドの媒介量と最終繁殖量を区別。

**Phase 3（分子／歴史論文）**：琉球白/淡紅と台湾白/青紫の花器官色素/発現＋ゲノム系譜＋同一背景での活性化・遺伝性を段階的に評価。反復的な表現型から「失活経路の再利用」を推定するのではなく、機能を測る。**時間を要する交配世代は独立の将来実験**とし、今すぐの結論にしない。

**Phase 4（一般化／統合）**：Q1 の別Cirsium亜節と北米系統のホールドアウト、そして**独立した別科（例えばchunのCamellia、Ho & SmithのIochrominae等）**で「形質×系統に依存するか」を見直す。3問はひとつの統一研究テーマになるが、査読向けには**機能論文・進化歴史/遺伝論文・他系統比較論文**を異なる主評価量に分ける方が無理がない。

## 「Nature級」の採否と反証例

一番強い、しかしまだ全く確認されていない主張候補：

> *An ancestral resource inventory can persist after a reproductive phenotype disappears; its later reassembly depends both on lineage-specific genetic accessibility and on which mutualists and antagonists confer fitness advantages to the resulting capitulum.*

**これが既知理論の言い換えで終わらない条件**：同一祖先的資源に由来する機能を分子的に確認し、異なる昆虫・環境条件での同一頭花モジュール組合せの適応度を測り、共通園・遺伝系譜・保留系統の将来反応で予測すること。ひとつでも未確認ならそのリンクは示唆として扱う。

**重要な反証**：同じ形質の記憶がCirsium内の別系統へ予測移植できても別科ではできない／真の刺は昆虫の接触を妨げず匂いが誘引を説明する／送粉と食害が同じ最適な頭花角度を好む／白花の色素構造遺伝子が失活していて実験的再発現できない／着色は浸透交雑か移入で、祖先の再発現ではない／資源は残るが遺伝的・発生的互換性がなく三資源を組み合わせられない。**いずれも価値のある結果**で、正方向のシナリオを出すためにゲートを緩めない。

## 2026-10-19 Q1の分子ソース修正：MN871982は対象種のITSではない

**配列番号・種名・ゲノム区画の対応を一次文献で精査**したところ、`MN871982` を `C. schantarense` とする後続2022年中国文献の引用は、元の[Xu, Dang & Li（2020, DOI 10.1080/23802359.2020.1721032）](https://doi.org/10.1080/23802359.2020.1721032)が明示する **別種 *C. shansiense*、甘粛省標本 WNU52012、葉緑体全ゲノム**と食い違った。この配列を対象 `C. schantarense` の核ITSとして使用しない。両種はKew POWOでも別の受容名である。

韓国の[Bae（2015）](https://doi.org/10.5352/JLS.2015.25.2.243)が比較した `C. shantarense/schantarense` 核18S–ITS1–5.8S–ITS2–部分28Sの**研究そのものは実在**するが、標本voucherとの同定を満たすアクセッションIDは現時点で**未回収**。既存の Moreyra2025多座位核標本2種と、2015単一座位資料を遺伝距離と見なして貼り合わせない。

**Q3にも直接関係する誤判定リスク：** [`C. schantarense f. albiflorum`](https://powo.science.kew.org/taxon/994340-1) という白花型の植物名は現在 Kew では `C. schantarense` の同物異名とされる。花色が白いことは、アントシアニン調節・構造遺伝子の共同保持、誘導可能性、表現型再獲得、遺伝性の実証ではない。

- [一次資料と二次引用の修正表](CIRSIUM_SCHANTARENSE_VS_SHANSIENSE_NUCLEAR_PLASTOME_ACCESSION_AUDIT_20261019.md)
- `analysis/validate_cirsium_schantarense_accession_provenance_v1.py` は種名＋voucher＋葉緑体／核マーカーの矛盾を自動拒否する。
- `analysis/probe_cirsium_schantarense_ncbi_nuccore_v1.py` は公共NCBI索引を読み取り、検索失敗を「DNA不存在」と呼ばず、NCBI候補を論文標本と同定したことにしない。

## 3問を単なる別々の論文にしない「2つの橋渡し検定」

**Bridge A：進化的記憶は何を覚えているのか（Q1×Q3）**

同じ白色または同じ頭花姿勢の表現型 `Z` を持つ複数のアザミ系統を比較し、異なる進化史 `H` が、独立検証された共同機能資源 `G` の有無および次の標準化環境・許可された発現試験での `Z_future` をどれだけ予測するかを試験する。

事前に区別する3つの予測器は `Z_visible` のみ、`Z_visible + H`、`Z_visible + measured functional G + H`。評価は**独立個体/集団・未使用系統の将来の発現結果**に対する予測損失とする。もし可視表現型しか測らないと歴史が有用でも、共同資源 `G` を独立に十分な精度で測ると歴史の追加予測力が消えるなら、**「記憶」は不可解な先祖の力ではなく、隠れた現在の資源に媒介される**という anto v9 の実証候補になる。歴史効果が残れば**測定していない現状態や環境**もあり得るため、非Markov性を証明したと呼ばない。

**Bridge B：利益があっても、なぜ形質が戻らないのか（Q2×Q3）**

相互作用相手の構成 `P` を変えると、ある頭花形態 `Z_1` に対して`Z_0`より大きな実測 `W` の利得が生じても、`Z_1` を形成する functional `G` が失われた系統にはその進化応答を期待できない。一方、`G` が保持されている系統では、適切な発現条件と遺伝伝達が確認できれば再出現する可能性がある。

これを検証するには**生態的価値（同一背景での `W(Z_1,P)-W(Z_0,P)`）と、遺伝的実現可能性（保有資源・発現・次世代伝達）を独立に操作/測定**する必要がある。系統間で `G` が異なれば、形質の適応度だけから「戻る速さ」を推定できない。逆に `G` が存在しても、互換性/発現が足りなければ戻らず、再出現しても昆虫相手の現在の環境で利得がなければ保持されない。

**統合候補の最強の反証可能な観測は：**「同じ形質状態の2系統が異なる機能資源を保ち、同じ相手群集変化に対して異なる発現可能性/子孫伝達を示す。同じ機能資源でも相手群集を変えると繁殖上の利得の方向が変わる」。これは現状どのリポにも**現実のCirsium観測としては存在しない**。従って大統一の先取りを禁止し、3問を別々に棄却可能にした上でBridge A/Bを**追加事前登録の統合検定**とする。

## 2026-10-19 source-backed decision update: Q1 actual external taxa, Q2 JPN36 conservation, Q3 leaf-only molecular states

The source-pinned audits now make the *actual* support for all three questions explicit, instead of treating planned sampling or another repo's simulation as original Cirsium biology:

| Question | Exactly observed/frozen now | What is still required | Current decision |
|---|---|---|---|
| **Q1 memory portability** | Japanese discovery roster **38**; NMNS states in **24**, definite head orientation **19** (13 up / 6 down), phyllary **10**, stickiness **13**. Independently sourced East Asian screen **19 taxa**, **16** non-overlapping; exact prequalified orientation **11** (6 up / 5 down), phyllary **1**, stickiness **2**. A taxonomically entangled `C. japonicum` and other discovery overlaps/ambiguities excluded. | Same source-validated nuclear phylogeny/branch distances including the true heldout taxa, enough replicated clades and homologous orientation values. Published 2-tip Arenicola pair cannot support a distance-dissimilarity correlation. | **Observed categorical trait diversity exists in both groups; Q1 memory portability remains NOT_IDENTIFIABLE.** Details: [real source audit](CIRSIUM_Q1_REAL_HELDOUT_TRAIT_COVERAGE_AND_TREE_GAP_V1.md). |
| **Q2 partner-dependent physical benefit** | EAzami JPN36 appressed phyllary and nonsticky history has the strongest localized candidate; original proposed 24-individual paired manipulation vs sham, **30 separate observation-only reference plants**. But **0 eligible/assigned/observed field heads** and no permission. Source cap of at most **2% flowering plants** implies **at least 1,200 counted flowering plants** to use all 24 experimental heads, subject to all separate permissions. | Document independent insect species/guild, real touch/contact/failed entry, sham integrity, confirmed egg/pollen and filled achenes in a real head cohort. | **Feasibility conditional and permissions absent, not evidence of no structural filtering.** |
| **Q3 reachability after disappearance** | The 2026 PRJNA1311153 registered metadata include **25 distinct real young-leaf RNA libraries**, among them 3 `C. brevicaule`, 3 `C. irumtiense`, and 3 white/3 blue-purple `takaoense`. All **18** predeclared 6×3 Ryukyu ANS/DFR/GMYC1 screening cells have **NOT_RETRIEVED/NOT_TESTED** functional candidate state. Phase-A **298 plants**, including 40 `C. sieboldii`, are planned **not enrolled**. | Same-individual flower-specific pigment chemistry, functional enzyme/regulated pathways, allele genealogies/introgression tests, inducible expression and heritability. | **No Cirsium regain, genomic LEGO mechanism, or fitness effect established.** |

**Crucial design distinction:** JPN36's proposed field 24 plants and Phase-A `C. sieboldii` 40 genetic plants are *not yet demonstrably the same plants*; do not add them as independent evolutionary histories or assume one includes the other. The current 11 outside-Japan orientation candidates have genuine non-overlapping **botanical state observations**, not automatically independent evolutionary origins. All three questions retain individual falsifiers and no new field authorization is implied.

**Next scientific action priority:** qualify a single genuinely shared head-and-plant ID, independent botanical orientation measurement and voucher-backed ancestry *before* crossing field access/fitness, and simultaneously seek a comparable multilocus placement for outside-Japan heldout orientation taxa. This is an empirical support problem, not a reason to invent additional simulated trade-off parameters.

## 2026-10-19 Q1追加の核DNA照合：最初の本物の予測対象はどこまであるか

権威的な形態データの**日本19分類群（上13／下6）対外部11種（上6／下5）**だけでは、系統的記憶の強さの移植は検定できない。新たにMoreyra 2025の**258行**の公刊Cirsium標本台帳を正確なSHAで照合したところ、外部11種のうち樹木コードとして同名が存在したのは**3種**、SRA提出名・BioSample/Runまで厳密に一致するのは**上向き2種（`C. vlassovianum`、`C. argyracanthum`）のみ**。下向きの `C. lidjiangense` という樹木コードにはトルコ産・提出名不一致の標本が該当し、中国の同名記載と安易にリンクできない。下向き0種の状態対照で `rho` を作ることは不可能。

これは**サンプルが存在するか否かについてのソース検証結果**であり、上向き/下向き状態のDNA配列・枝長・独立進化の再現性を推定したものではない。よって現状のQ1の最も価値ある次工程は、既存の下向き種の核ITSなどと同一のターゲット領域・樹木を確保できるかの精査である。ヘテロな遺伝マーカーの枝長を直接混合しない。詳細：`docs/CIRSIUM_Q1_REAL_HELDOUT_TRAIT_COVERAGE_AND_TREE_GAP_V1.md` と `data/evidence/cirsium_q1_real_japan_east_asia_holdout_trait_coverage_v1.json`。

## References and source ceiling

- [Chang et al. 2026, BMC Plant Biology 26:545, DOI 10.1186/s12870-026-08097-6](https://doi.org/10.1186/s12870-026-08097-6) (37 accessions, independent Ryukyu clades, Taiwan radiation, genome-size contrasts).
- [Ho & Smith 2016, DOI 10.1186/s12862-016-0675-3](https://doi.org/10.1186/s12862-016-0675-3) (non-Cirsium molecular intact anthocyanin machinery in unpigmented lineages; **possible regain** already documented, thus "gene retention makes regain possible" alone is not novel).
- [Larter, Dunbar-Wallis, Berardi & Smith 2018, DOI 10.1093/molbev/msy117](https://doi.org/10.1093/molbev/msy117) (Iochrominae repeated anthocyanin-expression shifts; candidate independent outgroup relevance).
- Existing signed head phenotypes and PR #44 observational contracts remain preserved, not merged into a causal result.
- **As of now: zero new same-individual real Cirsium genotype×flower-organ functional-retention×access×filled-achene linked records.** Existing 25 young-leaf RNA accession records and 2020 other-species gene family candidates do not meet this gate.
