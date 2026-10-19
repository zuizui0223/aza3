# *Cirsium schantarense* ITS 回収の障壁：MN871982 は別種の葉緑体ゲノム

**2026-10-19 | aza3 PR #50 | source-level correction to Q1; Q3 white-form pitfall.**

## 結論（実配列を取得したという主張ではない）

下向き頭花の *Cirsium schantarense* の核ITSは、**Bae (2015)** に原著実験として存在する。ただし2015年の韓国個体について、照合済みGenBank accession・voucher番号・解析可能なFASTAを今回の取得経路では確定できなかった。

核サンプル不足を埋めたように見える accession **MN871982** は、**Xu, Dang & Li (2020), DOI [10.1080/23802359.2020.1721032](https://doi.org/10.1080/23802359.2020.1721032)** により報告された **別種 *Cirsium shansiense*** の**葉緑体全ゲノム**（約152,396 bp）であり、採集地は中国甘粛省 Qingyang、標本 WNU52012。核ITSでも、韓国の *C. schantarense* でもない。

ところが[中国の別論文（2022）](https://html.rhhz.net/YXXB/html/20220735.htm)の本文には、`Cirsium schantarense` の葉緑体参照として **同じ MN871982** が記載されている。**一次資料と後続資料で accession と種名の対応が矛盾している**ことが、今回原著を突き合わせて判明した。後続の種名表記を正として元配列を自動的に付け替えてはならない。

Kew POWOは *C. shansiense* と *C. schantarense* を別々の**accepted species**として扱い、*C. schantarense f. albiflorum* は後者の synonym とする。白花型という学名が存在しても、それだけで「白花形質が消えた後に進化的経路が戻った」「独立した新系統が起源した」証拠にはならない。

## 重要な3階層を分離する

| 証拠 | 標本・由来 | ゲノム区画 | 結論上限 |
|---|---|---|---|
| Bae, Young-Min 2015: [10.5352/JLS.2015.25.2.243](https://doi.org/10.5352/JLS.2015.25.2.243) | 韓国 Dooryoon Mountainの *C. shantarense/schantarense* 1個体。*C. japonicum* 8地点などと比較 | 核rDNA 18S–ITS1–5.8S–ITS2–部分28S | 下向き種が上向き *C. japonicum* とITS系統で近い。**既知の相同配列研究**であるが accession 未回収、単一の連鎖rDNA領域 |
| Xu, Dang & Li 2020: [10.1080/23802359.2020.1721032](https://doi.org/10.1080/23802359.2020.1721032) | *C. shansiense* 甘粛省、voucher **WNU52012** | **plastid** 152,396 bp, accession **MN871982** | 全く別種の葉緑体。核系統の証拠に不可 |
| 2022年中国の[葉緑体研究](https://html.rhhz.net/YXXB/html/20220735.htm) | `C. schantarense, MN871982` と記載 | 引用された葉緑体参照 | **二次引用内の種名矛盾**。元論文の標本起源を上書き不可 |
| Kew POWO [*C. shansiense*](https://powo.science.kew.org/taxon/urn%3Alsid%3Aipni.org%3Anames%3A195998-1) と [*C. schantarense*](https://powo.science.kew.org/taxon/urn%3Alsid%3Aipni.org%3Anames%3A195961-1) | 同一種名ではない | 分類資料 | 文字列類似性、歴史上の綴り違い、var./f.を使って配列や祖先状態を付け替えない |
| Kew POWO [*C. schantarense f. albiflorum*](https://powo.science.kew.org/taxon/994340-1) | 白花型名の分類上のsynonym（原記載は1996年、タイプなし） | 色彩形態の**分類名情報** | 分子経路の保持・再活性化・遺伝性・再進化は未検証 |

**GenBank MN871982 の記録自体を、この監査で独立にFASTAとして取得・検証したと主張しない。** 一次論文が明示した accession と標本から検証した出典対応であり、配列照合は別ゲート。

## 今のQ1、Q3への具体的な帰結

**Q1：系統間記憶の転移。** 2015核ITSの局所系統近接性は実際の原著知見。ただし現在の公刊 Moreyra 2025 **同一多数核座位枠組み**への核標本照合は *C. schantarense* について0件（別の *C. vlassovianum* なら安全な標本照合がある）。2015の枝長を2025の枝長と混ぜて距離対非類似度 `rho` を作成しない。*C. schantarense* の2015核ITS **実際のGenBank登録番号を原著の全体・GenBank登録で別途確認**し、その標本と一致する多座位核系統のvoucher対応を探す。

**Q3：消失後の可能性。** 白花の `f. albiflorum` という名前だけでは、種内白花の頻度も、変異体の遺伝資源の可逆的発現も確定しない。この資料は白花表現型の探索候補にできるが、着色系統との家系、遺伝子型、小花の色素量・転写・色素再誘導を測る必要がある。葉緑体配列から核アントシアニン調節子の残存を推定しない。

## 今回の修正を機械的に検証する

- `data/evidence/cirsium_schantarense_vs_shansiense_primary_accession_audit_v1.json` — 出典ごとの種名・voucher・核/葉緑体区画・一次資料優先を明示。
- `analysis/validate_cirsium_schantarense_accession_provenance_v1.py` — *schantarense/shantarense* の綴り同定は対応可能、*shansiense* には**拒否**、核ITSとplastidには**拒否**、MN871982の用途は**原著が確定する *shansiense* のplastome以外に拒否**、synonymの白花学名を独立した再進化イベントに変換する場合も拒否。
- `.github/workflows/audit-cirsium-schantarense-accession-v1.yml` — 科学的誤結合の合成汚染を拒否し、**実配列未回収**、**遺伝記憶未計算**を保持する。

現実のDNA配列入手・種内分化・形態測定・進化歴史の新しい実証は **0**。これは**既存公開原著からの否定的な出典監査結果**であり、付加的な生物学的発見を捏造しない。
