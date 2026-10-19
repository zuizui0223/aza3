# 3問・第1問の実データ陰性対照：頭花状態の移植と進化的記憶の移植は別

2026-10-19 | `aza3` PR #50 | Q1 source-backed categorical orientation comparison, **NOT a new phylogenetic signal, rate or evolutionary reversion estimate**.

## 今回初めて直接照合した原本

Frozen **EAzami** commit `c9870e91604ed1514344ddc8f765eae5237dbf27` にある4点を、species/taxon labels を取りこぼさず入力した。

1. `data/evidence/moreyra2025_japan_38_membership_audit_2026-08-10.csv` — 日本に関係する38の公刊核系統tree内の taxon concept。
2. `data/evidence/japan38_nmns_capitulum_trait_seed_v1.csv` — 国立科学博物館（NMNS）の頭花状態22分類群。
3. `data/evidence/japan38_nmns_capitulum_trait_seed_extension_v2.csv` — 同ソースで追加した JPN34 と JPN24 の2分類群。上2つ合計24分類群で、`JPN` を38名簿へ厳密照合。
4. `data/evidence/heldout_east_asia_trait_state_audit_v1.csv` — Flora of China、2026年の韓国・中国分類記載などによる候補19種。そのうち日本の発見データと非重複で、事前定義された `candidate_primary` 頭花向きは11種。

本レイヤーは**独立した観察画像/配列を新規に測定したのではなく、権威ある分類形態記述の実在する状態ラベルを比較している**。候補11種は日本の探索分類群と重ならないが、**11回の独立進化起源でも、11個の系統学的に独立した反復でもない**。

## 実際のsource-native照合結果

| 判定 | 日本Japan38の公刊名簿・NMNS分類形態 | 東アジア独立候補 |
|---|---:|---:|
| 原本の全分類群 | 38 | 19候補 |
| 状態ソースがある | 24 | 非重複16 |
| **厳密に調和した頭花向き二状態** | **19** | **11** |
| 上向き／直立 | 13 | 6 |
| 下向き／垂下 | 6 | 5 |
| 総苞片姿勢の分類形態記録 | 10 | 1 |
| 粘着性の分類形態記録 | 13 | 2 |

日本のstrict向き19は `upward_or_erect` / `downward_or_nodding` の**完全同じ定義**に限る。元ソースの `JPN_32=upward_or_ascending` を自動的に「upward_or_erect」とみなさず、`JPN_17=source_conflict_index_downward_detail_erect` を二択へ強制しない。さらにunknown3を陰性・反対状態として扱わない。したがって元 `EAzami` が文脈によって報告した「向き約20の解決taxa」と、この**厳密共通コード19**は矛盾しない。

東アジア19候補からは、探索に重複するもの、既存名と分類学的に絡むもの、`ambiguous` の4行、未知／同定不能状態を**先に定義されたステータスで除外**する。11種の向きは上向き6（例：*C. vlassovianum*, *C. vernonioides*）、垂下5（例：*C. reflexiphyllarium*, *C. fangii*）。**自然集団個体内の定量角度ではなく、分類記載でのカテゴリー**である。

外部の `nonsticky` と日本の `nonsticky_or_nearly_nonsticky` は測定語義が完全に一致せず、`recurved` と `spreading_or_recurved` も同じコードへ自動変換しない。ゆえに今回**比較可能なゼロショット二状態予測は頭花向きだけ**で、総苞・粘着では類似状態の「科・亜節間で記憶が保存された」という結論を出さない。

## まだ系統を使わない、明示的な負の対照

日本のsource19から推定される頭花が上向きの比率は

`p_J = 13/19 = 0.6842105263`。

この一定確率だけで東アジアholdoutの `y=1`（上向き6種）、`y=0`（垂下5種）を予測すると：

- 未使用東アジアの実際の上向き比率は `6/11=0.5454545`。
- 日本の優勢状態「上向き」を全外部11種に割り当てると **6/11、54.5%** で一致。
- 平均Brier（二値確率の二乗誤差） `mean[(y-p_J)^2] = 0.2671871065`。
- 何も知らない等確率 `p=1/2` のBrierは `0.25`。
- `0.2672 > 0.2500`：**このsource sampleと単一の周辺頻度予測では改善しない**。

これは有意差検定でも、独立進化イベントのn=11でもない。東アジアのラベルをすべて見た後で求めた `p_{EA}=6/11` を校正比較に使うことも、予測性能の評価としては不適切。分類記述の発表/選抜、種名の公刊歴と環境が未調整であり、これを「アザミの進化的記憶は系統間に通じない」という結論に**絶対に格上げしない**。

**むしろ今回確認できた肯定的な事実**：日本の発見パネルと非重複な東アジア候補で、頭花の上/下状態の両方を独立の植物誌・種記載ソースから回収できた。これは真にout-of-sampleな**頭花状態ラベル照合の入り口**であり、進化記憶の外部検証を実装するための必要な材料の一部である。

## 公刊核系統の標本記録まで照合：11種のうち本当に接続可能なのは2種

さらに `EAzami@c9870e...` に保存された **Moreyra et al. 2025 のSupplementary Table S1核遺伝データ標本記録258行**（8分割CSV＋原本SHA256 manifest）と、前述の東アジアの向きが確定した11種を**種名、BioSample、voucher、地域、樹木コード対SRA提出名の整合性**で照合した。

| 外部候補 | 公刊核サンプル記録 | 判定理由 | 頭花の向き |
|---|---|---|---|
| *Cirsium argyracanthum* | **SAMN44017818**, チベット・中国、voucher PE01890882 | S1 tree codeと提出名が一致、原標本は中国、`core_east_asia` | 上向き |
| *Cirsium vlassovianum* | **SAMN34240350**, モンゴル、voucher MW0194150 | S1 tree codeと提出名が一致、原標本はモンゴル、`northeast_asia_bridge` | 上向き |
| *Cirsium lidjiangense* | SAMN34240312 | **危険な種名だけの一致**。原標本はトルコ・Bursa、SRA提出名とも食い違い、地域区分も `outside_scope`。**外部DNA支持として不採用** | 植物誌では下向き |
| その他の外部候補 | 当該S1標本表で該当tree codeなし：**8種** | 未回収はDNA不存在の証明ではない。別公刊系統/同物異名の照合は別途必要 | 上下混在 |

したがって、今回のsourceで確認した数は **可信な核サンプル一致2種／名前のみの危険な一致1種／一致なし8種**。この2種は同じグローバル核研究の**標本登録**に接続できるため、「東アジア11種に核データの接続が一つもない」という扱いは不正確だった。反面、この研究の完全な最終系統樹Newickと比較可能な枝長は、このソースでは回収できていない。**標本表での一致＝そのサンプルを同じ系統樹へ置いてpatristic距離を計測済み、ではない。**

**さらに決定的な問題：信頼できた外部2種は両方とも上向き**。外部2種だけの唯一の形質ペアは非類似度0であり、たとえその対の枝長があっても、外部系統について系統距離と形質非類似度の勾配 `rho` を推定できない。下向きの外部候補の核樹上位置と出典に一致した状態を少なくとも含める必要があるが、必要な系統分枝の独立反復は別途検出力ゲートを設定する。

ここは新しい**実データに基づく到達可能性の診断**である。改めて、これは遺伝的部品の保持、トゲ・向きの適応、trait-history portabilityの肯定的な測定結果ではない。

- 実行コード：`analysis/audit_cirsium_q1_heldout_nuclear_sample_crosswalk_v1.py`
- 厳密な原本：`EAzami/data/evidence/moreyra2025_cirsium_reconciliation_v1/manifest.json`、原サンプル258行＋8部分ファイルのSHA256。
- CI：`.github/workflows/audit-cirsium-q1-nuclear-heldout-sample-v1.yml` が同じソース版を取得して、`C. lidjiangense` のトルコ/名称矛盾の見落としを拒否する。
- **次の実行優先順位**：11候補を新たに配列取得する前に、他の公刊核データに保存された標本コード、tree tip、voucher、SRA提出名を検索し、**下向き候補の実際に比較可能な枝長付きの系統標本**を探す。上向き候補の数を増やすだけでは外部記憶勾配の識別性は上がらない。

## 2026-10-19 分布点品質 × 核標本 × 向きの三重交差

本研究で必要なDNAデータの欠損位置は、元の `EAzami` [strict occurrence QC](https://github.com/zuizui0223/EAzami/blob/c9870e91604ed1514344ddc8f765eae5237dbf27/docs/chapter2/HELDOUT_EAST_ASIA_OCCURRENCE_SUPPORT_RESULT_V1.md) と照合すると、さらに明確になる。

| 独立候補種 | 公刊頭花向き | 条件を通ったGBIF分布点 | 同じMoreyra核研究の安全なBioSample |
|---|---|---:|---|
| *C. schantarense* | 下向き・垂下 | **8** | **なし** |
| *C. vlassovianum* | 上向き・直立 | **8** | **あり：SAMN34240350** |
| *C. argyracanthum* | 上向き・直立 | **0**（厳格な座標精度条件） | あり：SAMN44017818 |
| その他8候補 | 上・下が混在 | 主データとしての条件に不足 | 安全な完全一致なし |

GBIF監査は「正確な分類名・座標誤差10 km以下・0.1度間引き・少なくとも3点」を対象とした、**先に固定された出現データの品質試験**である。環境ラスタのBIO値は開かれていない。分布点8件は、独立な系統変化8回ではない。

よって外部11種中、
- 明瞭な向き + 安全な核サンプル：**2種**（ともに上向き）。
- 明瞭な向き + 質の高い分布点：**2種**（上1・下1）。
- 三者が揃う：**1種、上向きの *C. vlassovianum* のみ**。
- 三者が揃う下向き種：**0種**。

**次の具体的な遺伝学上のボトルネックは、下向きで分布点が8件ある *C. schantarense* の、既存の公刊核標本の正確な同物異名／voucher照合、あるいは同一遺伝子座群に接続できる核配列を確認すること**。ただし1種を追加しても進化的記憶の `rho` を推定する独立な系統反復にはまだ不足するので、そこで統計量や自然選択を先取りしない。

原本は `EAzami@c9870e...` の `heldout_east_asia_occurrence_support_v1.csv`（12分類群のスクリーン）であり、新設 `analysis/audit_cirsium_q1_heldout_nuclear_sample_crosswalk_v1.py` にGBIF品質との交差監査を追加した。小さい[機械可読の確定レシート](../data/evidence/cirsium_q1_nuclear_occurrence_source_overlap_v1.json)に、2/11、2/11、1/11と不採用になった危険なトルコ名一致を別々に保存した。

## 重要な補足：下向き *C. schantarense* の核DNA研究自体は既に存在する

「Moreyra et al. 2025 の**共通Compositae1061標本表に下向き候補が整合しない**」ことと、「*C. schantarense* の核DNAがまったく研究されていない」ことは別である。

**Bae（2015）, *Genetic Relationship of Some Cirsium Plants of Korea*, *Journal of Life Science* 25:243–248, DOI [10.5352/JLS.2015.25.2.243](https://doi.org/10.5352/JLS.2015.25.2.243)** は、韓国の垂下した `C. schantarense`（原文中に `shantarense` の綴りも見られる）1個体と、上向き *C. japonicum* 8地点などについて、**18S–ITS1–5.8S–ITS2–一部28S の核rDNA**を比較した。原著の隣接結合法のITS系統では、向きが逆なのに両者は近い群に位置した。これは既に `EAzami` の独立した[文献Tier C監査](https://github.com/zuizui0223/EAzami/blob/c9870e91604ed1514344ddc8f765eae5237dbf27/docs/CIRSIUM_PHYLOGENY_LITERATURE_SCREENING_BATCH_02B_2026-08-10.md) に収録されていた結果であり、ここで新発見を主張するわけではない。

**非自明な生物学的な示唆候補：** 目立つ頭花の向きの差は、浅い核リボソームDNA変異でも起こり得る。しかし `ITS` は**一つの連鎖した反復遺伝子領域**で、concerted evolution・浸透交雑・ILS・微小な遺伝差を区別できない。垂下頭花が最近再進化した、祖先遺伝子が復活した、形質の記憶が消えた、という証明にはならない。

- 2015年の配列原本・GenBank番号・当時の標本voucherを**今回の解析では再取得していない**。
- 2015年のlinked rDNA系統と、2025年の核1061遺伝子座系統について、異なる枝長を貼り合わせて `rho` を算出することは禁止。
- 次の解析は *C. schantarense* の**DNAが存在するかどうか**を初歩から尋ね直すのではなく、すでにある2015年の核ITS資料の標本同定・配列accessionを監査し、可能なら2025年の多遺伝子座共通核枠組みにどのように接続できるかを問う。

## 進化的記憶を検証する追加条件

元の `turnover` の量は、系統のpatristic距離と同一形質の非類似度の関係（S3 rhoなど）。**今回の `p_J`, Brierにはpatristic距離が1つもない。** 正確なQ1は次を満たす場合に限る。

- 日本38と外部候補11を含み、taxon concept・voucher・核遺伝子樹上の位置が一致した**同じトポロジー/枝長尺度**、ならびに topology ensemble（交雑・細胞型を考慮）。
- 少数の極端に近縁な組み合わせに支配されない、十分な独立な部分系統とtrait状態。たとえば琉球 `Arenicola` の2種だけでは一対しかなく、記憶 `rho` を推定できない。
- 同じ発育段階で相同性ある、調和済みの頭花向きと総苞・粘着性の現物測定（写真の画像上下は重力方向ではない）。
- 形質・地域・分類群を**見てから樹木を都合よく選ばない**。日本で学習、外部系統を丸ごと保留し、距離・特徴の予測改善を同じ損失関数で比較。ペア数 `n(n-1)/2` を独立サンプル数と誤認しない。
- 属内で反復できても、**同じ形質の記憶が別の植物科に通じる**には、さらに相同形質の別科データ・独立系統の比較が必要。

このうち核系統距離の外部候補11種の一貫したマッチングがないため、現在の正式な Q1 判定は `NOT_IDENTIFIABLE_MISSING_MATCHED_PHYLOGENETIC_DISTANCE_ENSEMBLE`。

## 再現コードと厳密な結果の地位

- `analysis/audit_cirsium_three_question_trait_state_portability_v1.py`：sourceの38→24 join、解決不十分なJPN17/JPN32、source-conflicted state排除、東アジア19→非重複16→向き11、カテゴリー分布、Brier計算、系統記憶HOLDを出力。
- `.github/workflows/audit-cirsium-head-state-transfer-v1.yml`：毎回EAzami commit `c9870e...` の**実物4原本**を取得し、失敗状態7件（重複、未確認状態を上向きへ昇格、系統名錯誤、発見/保留漏洩など）をsynthetic corruptionで検証し、元データを変更しない。
- 状態頻度の診断は**実データの記述解析**。これを超える進化的記憶の数値や新たな再進化回数、分子的祖先の再利用、Q2の相手別利得、Q3の再発現の結果はまだ**未測定・未実証**。

**研究上の次の真のデータ作業**：外部11分類群の樹木上の正確なvoucher/tipを、東アジアの公刊核系統データと1対1で照合する。外部候補の樹木情報が一部のみなら、全てを一つの既知の系統に置いたような樹木を作らず、計算可能な下位集合だけを報告する。これは source portabilty と phylogenetic memory portability を識別する最短のゲートとなる。
