# アザミ頭花：匂い・時刻・発育段階の選別と、実際の侵入口の選別を分ける

**2026-10-08 | PR #44 | prospective inference boundary, not a functional field result.**

## 公刊されている「すでに存在する選別」

- **Theis 2006**, *Journal of Chemical Ecology*, DOI [10.1007/s10886-006-9051-x](https://doi.org/10.1007/s10886-006-9051-x): *Cirsium arvense* で13種中10種の主要花香成分を餌付きトラップで試験。ベンズアルデヒドとフェニルアセトアルデヒドを含め、送粉候補と花食者双方を誘引する成分、比較的片側を誘引する成分を確認。これは**化学信号が送粉者だけに専用ではない**という直接実験であり、頭花の総苞・刺の物理障壁の効果を測ったものではない。
- **Theis, Lerdau & Raguso 2007**, *International Journal of Plant Sciences* 168:587–601, DOI [10.1086/513481](https://doi.org/10.1086/513481): *C. arvense* で日周花香放出は送粉者の活動時間と対応し、花食者の多い時間には低い。一方 *C. repandum* は香りのピークが午後の日中で昆虫活動と必ずしも対応しない。生殖成熟期には両群の来訪・香りが増加。**日周パターンは観察上の対応であり、花香の周期が選択によって進化した証明ではない**。
- **Vanbergen et al. 2007**, *Ecological Entomology*, DOI [10.1111/j.1365-2311.2007.00885.x](https://doi.org/10.1111/j.1365-2311.2007.00885.x): *C. palustre* 240個体・24ブロック移植で、発育段階が *Tephritis conura* と *Pteromalus elevatus* の個体数・寄生率の変動と結びつく。隔離距離よりフェノロジーが重要な対抗説明。

**Implication:** 頭花の外形・向きと種子食者の訪花が相関しても、同じ種の昆虫が同じ日に頭花へ来るかは、香りと昆虫の日周活動によって変わる。よって未調整の侵入経路の違いを**真のトゲの選択的防御**と解釈できない。「花の時刻別防御」一般も先行研究であり、今回の新規性にはできない。

## 重要なデータの欠落：既存のイベントには絶対時刻がない

現在の `capitulum_video_effort_denominator_v1.csv` と `capitulum_guild_access_event_ledger_v1.csv` は頭花×観察boutの**相対秒数**で記録し、`phenological_stage`はあるが、採録当日の日時・日周時刻がない。したがって既存 `stage_matched_equal_stratum_overlap` は、異なる日の朝の送粉者・夕方の食害者を「同じ開花段階の異なる入口」と誤認し得る。これは実際のアザミで確認された誤りではなく、解析上の選択・交絡リスクである。

既存頭花・boutの列は一切変更しない。新たに **空の** `data/intake/capitulum_clock_matched_effort_sidecar_v1.csv` で固定する：

- `individual_id, population_id, capitulum_id, observation_bout_id`: 既存の真正joinキー。
- `video_start_local_iso, video_end_local_iso`: タイムゾーン付き実時刻。1区間 **15分以下**。検証済みの頭花+接近領域の連続映像と絶対的に同じ時間長が必要。
- `recording_evidence_uri`: 実映像の検証参照。`contemporaneous_weather_reference`: 外部気象記録の参照または明示的 `not_assessed`。
- `floral_scent_status`: `measured` / `not_assessed`。前者は独立の `floral_scent_sample_reference` が必要。時刻が合うだけでは**匂いを統計的に調整できたとは言わない**。

## 実装：同時・同段階・同頭花順位でない限り比べない

`analysis/audit_capitulum_concurrent_clock_guilds_v1.py` は、既存のビデオ努力量検証（視野、event-trigger-only排除、真のゼロ、現場ID、ステージ）を通してから、`population × stage × head rank × reproductive sex state × exact calibrated UTC start/end` でマッチする。15分間の**同時刻連続撮影**であることが必要。もし一方のギルドがいなければそれは「比較可能な侵入経路差」のデータではない。

有効なギルド別結果には最低限、**入口前に独立に同定**した各10接近、各3独立植物、未分類経路0。これは計算可能性の最低条件で、検出力/植物集団の独立進化起源を示すものではない。観測映像から`reproductive_zone_reached`が分かれば、区間別の成功回数/独立接近数も記録するが、未評価`NA`を失敗0にしない。複数イベントは植物あたりで相関し、人工的にサンプル数を稼げない。

**さらに厳しい反証可能性：** 一つの開花段階で午前と午後に、送粉者・食害者が**各時間窓内では全く同じ世界方向の経路**を使っていても、午前90:10、午後10:90というギルド組成の違いにより、段階全体をプールすると重複度は `0.2`、完全に一致した同時刻を比べれば `1.0` となる。この例は**純粋な合成データであり、野外観測ではない**。解析のCIテストはこの偽の棲み分けと時刻不明/写真のみ/ギルド事後付与を拒否する。

## 因果的なメイン問いを一つに絞る

**Primary estimand:** 同一種・同程度の到着圧・同じ頭花発育段階/日周時間/生殖状態で、実測 `natural_orientation_deg × spine_length_mm × spine_rigidity × phyllary_gap_mm` の変異が、実際の**失敗を含めたアクセス成否**に及ぼす影響。種子食者と正当な送粉者と寄生蜂では個別に評価する。

各部位の解釈境界：

1. **誘引/到来:** 視覚提示・花香・別種の開花資源・天候による潜在的到着率。到来分母が不明な録画を真のゼロとして扱わない。
2. **幾何的アクセス:** 本物の総苞・トゲへの最初の接触、阻害、迂回、最終侵入。同一種の昆虫・同時刻/同段階で比較し、花香だけの説明に対抗させる。
3. **機能:** 柱頭への実際の花粉沈着/産卵確認/被害前の寄生蜂による有効宿主死亡。目レベルのMobileNetV3は候補発見まで。
4. **最終適応度:** 成熟した充実・発芽可能痩果と（可能なら）花粉輸出。行動経路の差が直ちに自然選択になるわけではない。

**一段先の競合予測:** 香りの日周パターンがギルドの到着を分けたが、到着後には植物の総苞構造を素通りするなら、形態的トレードオフは機能として存在しないかもしれない。逆に同時刻・同段階・同種の昆虫が到着した後に、ある形態で種子食者だけ失敗するなら、物理的選択性が候補になる。ただし向き自体が香りの拡散や着陸位置を変える可能性は残り、最終的な因果分離には許可された操作か強い自然実験が要る。

## Current evidence and STOP

- Existing **Azami** 1,734 image observations/42 taxon labels: partial phenotypic coupling, not measured physical spines/gravity orientation. Phenotype context interactions BH支持0/12。
- GloBI: partner discovery complete, but source/target rows do not identify same-head timing or access.
- Existing `aza3` real continuous-bout / clock-sidecar / verified guild entry: **zero rows**. The new exact-concurrent model prints `NO_REAL_HEAD_VIDEO_OR_CLOCK_EFFORT` with null effect and no p-values.
- **STOP** if no overlapping guild/time support, if pre-entry taxon cannot be verified, if there is no head-stage validation, or if functional endpoints remain unmeasured. Do not turn the missed data into an adaptive zero.
- **No field experiment, sampling, manipulation or camera recording newly authorized by this code change.** Existing EAzami/Azami contracts and frozen results are unchanged.

**Next evidence milestone:** from a permitted, clearly observed flowering head, identify exactly one candidate species/behavioral role, the initial contact surface, true anatomical gap, whether access failed, and linked video/time reference. Until then, no positive claim about adaptive architecture.
