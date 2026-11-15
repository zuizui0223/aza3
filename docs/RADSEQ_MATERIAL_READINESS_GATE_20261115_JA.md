# RAD-seq 発注前の試料・形質対応ゲート（2026-11-15）

**Status: HOLD — 自前試料・提供試料の登録と受領を確認できていない。**

**RAD中心への変更案（2026-11-15）:** [日本産アザミRAD-first系統推定設計](RAD_FIRST_JAPAN_CIRSIUM_PHYLOGENY_V1_20261115.md)。ここに記した旧GROWのtarget-capture主従関係は比較資料として残すが、新案では新規核DNA取得をRAD優先とする。試料受領・許可・倍数性・RAD共通座位の停止条件は引き続き有効。
この文書は試料採集、実験、シーケンス発注の許可ではない。

## 結論と優先順位

1. **現行 GROW D1 主線**: 頭花の向き O と総苞構造 B を中心とした系統史。既存 Moreyra 350核座位の系統を広域の参照に用い、RAD は認可済み近縁集団・個体の祖先構造と交雑診断に使う。まず受領・voucher・個体ID・権利・DNA品質・倍数性を確認する。
2. **予備試験は24–48試料の暫定範囲**。保存法、2C/倍数性、制限酵素と共通座位、dropout、paralog/homeolog、技術反復一致率を、形質結果を開く前に評価する。試料未受領の状態で「24試料選定済み」と書かない。
3. **広域種系統の不足**は Compositae1061互換 target capture / 既存核DNAを使い、広範なRAD単独樹を強制しない。倍数性を無視した共通diploid SNP matrixは作らない。
4. **同一個体の実測**を原則とし、重力基準の頭花角・総苞の開き／真正トゲ長・配置・粘着性・花色と、DNA・voucher・cytotypeを個体IDで紐づける。歴史的な分類群コードの形質を新規DNA個体の実測値に置換しない。

## 台帳の実監査（GitHubの公開ファイル、2026-11-15）

| 対象 | 実行時の公開台帳行数 | 意味 |
|---|---:|---|
| [aza3 Nature個体台帳](../data/intake/aza3_nature_individual_intake_v1.csv) | **0** | 予定個体の実試料登録はない |
| [aza3 Chapter 3 個体台帳](../data/intake/chapter3_individual_intake_v1.csv) | **0** | これは従来のハンドオフ用空台帳 |
| [東大提供予定試料の受領台帳](../data/intake/grow_ito_sample_handoff_v1.csv) | **0** | 提供予定であり受領済みではない |
| [EAzami DNA/表現型台帳](https://github.com/zuizui0223/EAzami/blob/main/sampling/aim13_individual_sample_ledger_v1.csv) | **0** | 別repoの従来採集計画の台帳 |

**重要:** 台帳上0件と研究室の物理試料0個体は同義ではない。未登録の現物試料が存在する可能性は残る。公開リポに正確な所在地、収集許可証の本文や個人情報を記載しない。現地の資料照合は非公開の保管台帳で行う。

## 異なる採集計画を混同しない

| レイヤー | 数 | status |
|---|---:|---|
| [EAzami旧Core population計画](https://github.com/zuizui0223/EAzami/blob/main/sampling/SAMPLING_DESIGN.md) | 14集団・190個体 | 先行計画。現行GROWの必須採集数ではない |
| [aza3 Phase A genotype×phenotype構想](../docs/AZA3_NATURE_GENOMIC_SAMPLING_PANEL_V1.md) | 6種・298個体 | 広い形質モジュール比較の構想。現行GROW D1のpilot数ではない |
| [GROW D1集中3年計画](GROW_PHD_FOCUSED_ORIENTATION_PHYLLARY_3YEAR_V1_20261027.md) | 24–48試料の予備試験 | 現行の試験範囲。受領した試料の質によって可否判断 |
| [aza3 PR #3 全国マスターv11案](https://github.com/zuizui0223/aza3/pull/3) | 314 unique source のactive minimum、451 prospective slots | **未マージの提案**。228 tree slotsは野外で228個体を破壊採集する指示ではない |

2つの中心形質O/Bに絞ったGROWと、全球〜全日本を拡張するaza3の研究構想の「採集数」は足し算しない。

## 購入・外注前の4ゲート

| ゲート | 何を照合するか | 未達の扱い |
|---|---|---|
| G0 受領／権利 | 現物、採集・移転・解析・配列・共有の許可、保全審査 | **HOLD** |
| G1 ID／形質 | DNA個体IDとvoucher、実測O/B／撮影像の同一個体性 | 個体レベル解析に入れない |
| G2 DNA／倍数性 | 組織保存、DNA量・分解、フローサイトメトリー、cytotype | pilotライブラリの投入を保留 |
| G3 RAD技術pilot | 24–48試料（あれば）、酵素・サイズ選択・共通座位・技術反復・allele dropout・paralog | 共通遺伝座/倍数性で失敗ならRADを集団内に限定 |

**GOは自動判定しない。** 予備試験は候補試料が24未満なら不足と報告し、24以上になっても実試料・許可・DNA・倍数性・多型・対照の代表性を責任者が確認して初めて発注可能とする。

## 機械可読の現状確認

```bash
python analysis/audit_radseq_pilot_intake_v1.py
# 必要な申請用 documentary cohort が24個体に満たない場合、終了コード2
python analysis/audit_radseq_pilot_intake_v1.py --require-documentary-pilot-floor
```

この監査はaza3の**3つの公開台帳**だけを読み、識別子・正確な地点・許可資料をstdoutへ出力しない。G0–G3の外部物理確認は代替しない。陽性数は「文書上のゲートを満たした引渡し個体」であり、ユーザー・受入側の実務確認なしには生物学的な「確保個体」と称さない。

## 次に新しいデータとして必要なもの

**最初の1件は東大提供試料について、機微情報を除いた受領・権利・個体/標本対応・DNA保存/品質の事実を記入すること。** 受領がなければ現物照合・予備実験は進めず、既存の核系統樹と真正形質の対応づけを先行する。

専門領域のgate:
- [GROW D1 受領判定](GROW_D1_SAMPLE_HANDOFF_READINESS_AND_ALTERNATIVES_20261026.md)
- [現行3年計画](GROW_PHD_FOCUSED_ORIENTATION_PHYLLARY_3YEAR_V1_20261027.md)
- [Chapter 3 主線 issue #11](https://github.com/zuizui0223/aza3/issues/11)
- [P02 JPN06–JPN15 RAD falsification issue #8](https://github.com/zuizui0223/aza3/issues/8)

**限界:** RADの集団構造や系統位置のみから原因遺伝子・向きの機能・トゲの防御適応・再発現・独立進化回数の確証は得られない。各形質の機能・繁殖適応度検証は別の因果レーンに残す。
