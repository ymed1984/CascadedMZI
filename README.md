# CascadedMZI

Mach-Zehnder interferometer (MZI) / asymmetric MZI (AMZI) を使った光フィルタの簡易設計・スペクトル計算ツールです。

CLI から構成とパラメータを指定し、設計レポート、透過スペクトル CSV、スペクトル PNG を出力できます。

## Requirements

- Python 3.12+
- uv

依存パッケージは `pyproject.toml` / `uv.lock` で管理しています。

```bash
uv sync
```

## Quick Start

ヘルプを表示します。

```bash
uv run python main.py --help
```

Simple AMZI モードを実行します。

```bash
uv run python main.py
```

CSV と PNG を出力します。

```bash
uv run python main.py \
  --architecture simple \
  --csv out/simple_amzi_spectrum.csv \
  --plot-file out/simple_amzi_spectrum.png
```

論文参照の LAN-WDM demux モードを実行します。

```bash
uv run python main.py \
  --architecture lan-wdm \
  --span 45 \
  --points 451 \
  --csv out/lan_wdm_paper.csv \
  --plot-file out/lan_wdm_paper.png
```

`pub_5335.pdf` の 16-channel cascaded MZI WDM 例を実行します。

```bash
uv run python main.py \
  --architecture pub-16wdm \
  --points 801 \
  --csv out/pub_16wdm.csv \
  --plot-file out/pub_16wdm.png
```

## Architectures

### `simple`

複数段の AMZI を直列に接続した簡易モデルです。

各段は以下の構成として扱います。

```text
input -> coupler -> arm delay dL -> coupler -> output
```

各段の FSR は次のように決まります。

```text
stage 1: base_fsr
stage 2: base_fsr * fsr_scale
stage 3: base_fsr * fsr_scale^2
...
```

出力は 1 本の透過スペクトルです。各段の AMZI 透過を掛け合わせます。

### `lattice`

カプラと遅延線を交互に接続した MZI lattice filter モデルです。

```text
coupler -> delay -> coupler -> delay -> ... -> coupler
```

2x2 カプラ行列を電場ベクトルに順番に掛け、through port と cross port のスペクトルを計算します。

利用できるプリセット:

```text
dwivedi-flat
maxflat-n2
maxflat-n3
maxflat-n4
cheby-n4
```

例:

```bash
uv run python main.py \
  --architecture lattice \
  --lattice-preset maxflat-n4 \
  --csv out/lattice.csv \
  --plot-file out/lattice.png
```

### `lan-wdm`

`info/photonics-09-00252.pdf` を参照した、3-stage cascaded MZI 型の 8-channel LAN-WDM demux モデルです。

論文の Figure 1 の binary-tree 構成と Table 2 の遅延長設定を元にしています。

デフォルト値:

```text
center wavelength : 1291.0 nm
effective index   : 2.93
group index       : 3.83
channel spacing   : 4.5 nm
coupler ratio     : 0.50
```

ターゲットチャンネル:

```text
ch1 : 1273.5 nm
ch2 : 1277.9 nm
ch3 : 1282.3 nm
ch4 : 1286.7 nm
ch6 : 1295.6 nm
ch7 : 1300.1 nm
ch8 : 1304.6 nm
ch9 : 1309.1 nm
```

ポート対応:

```text
Port 1 : ch1 = 1273.5 nm
Port 2 : ch9 = 1309.1 nm
Port 3 : ch3 = 1282.3 nm
Port 4 : ch7 = 1300.1 nm
Port 5 : ch2 = 1277.9 nm
Port 6 : ch6 = 1295.6 nm
Port 7 : ch4 = 1286.7 nm
Port 8 : ch8 = 1304.6 nm
```

### `pub-16wdm`

`info/pub_5335.pdf` の Fig. 5 / Fig. 6 を元にした、16-channel cascaded MZI WDM demux の暫定モデルです。引用 [13] の `info/oe-21-10-11652.pdf` に記載された binary-tree CMZ の考え方に従い、段が進むごとに MZI のアーム長差を半分にします。

構成は 4 段です。

```text
MZI 1
  -> MZI 2.1 / MZI 2.2
    -> MZI 3.1 ... MZI 3.4
      -> MZI 4.1 ... MZI 4.8
```

デフォルト値:

```text
center wavelength       : 1550.0 nm
effective index         : 2.40
group index             : 4.20
channel spacing         : 800 GHz
coupler ratio           : 0.50
coupler excess loss     : 0.03 dB/coupler
default wavelength span : 110 nm
```

ポート対応は `pub_5335.pdf` Fig. 5 に合わせています。

```text
Port  1 : ch1
Port  2 : ch9
Port  3 : ch13
Port  4 : ch5
Port  5 : ch3
Port  6 : ch11
Port  7 : ch15
Port  8 : ch7
Port  9 : ch2
Port 10 : ch10
Port 11 : ch14
Port 12 : ch6
Port 13 : ch4
Port 14 : ch12
Port 15 : ch16
Port 16 : ch8
```

将来 Lumerical のカプラ設計結果を使う場合は、CSV を `--coupler-file` で指定できます。デフォルトでは以下の列名を読みます。

```text
wavelength_nm,coupling_ratio,excess_loss_db
```

例:

```bash
uv run python main.py \
  --architecture pub-16wdm \
  --coupler-file data/bent_dc_lumerical.csv \
  --csv out/pub_16wdm_lumerical.csv \
  --plot-file out/pub_16wdm_lumerical.png
```

列名が異なる場合は、以下で指定できます。

```bash
uv run python main.py \
  --architecture pub-16wdm \
  --coupler-file data/bent_dc_lumerical.csv \
  --coupler-wavelength-column wavelength_nm \
  --coupler-ratio-column coupling_ratio \
  --coupler-loss-column excess_loss_db
```

現在はカプラの暫定値としてフラットな 50:50 と `0.03 dB/coupler` を使っています。これは `pub_5335.pdf` の bent DC の説明に合わせた仮値です。実際の bent/straight DC 比較を行う場合は、Lumerical や測定から得た波長依存の coupling ratio を CSV として入力してください。

位相オフセットの効果を確認したい場合は、`--no-phase-offset` で自動位相合わせを無効化できます。

```bash
uv run python main.py \
  --architecture pub-16wdm \
  --no-phase-offset \
  --csv out/pub_16wdm_phase_offset_off.csv \
  --plot-file out/pub_16wdm_phase_offset_off.png
```

## Main Equations

FSR から遅延長を計算します。

```text
dL = lambda0^2 / (ng * FSR)
```

`pub-16wdm` のように周波数間隔で設計する場合は、以下を使います。

```text
dL = c / (ng * FSR_frequency)
```

論文参照の `lan-wdm` では、1段目の FSR をチャンネル間隔の 2 倍として扱います。

```text
dL1 = lambda0^2 / (2 * channel_spacing * ng)
```

チャンネル位置合わせ用の追加遅延は以下です。

```text
dLshift = lambda0 / neff
```

AMZI の位相は、中心波長まわりの一次近似として以下の形で扱います。

```text
phase(lambda) =
  phase_offset
  + 2*pi*ng*dL*(1/lambda - 1/lambda0)
```

## About `phase_offset_rad`

`lan-wdm` モードでは、各 MZI 段に `phase_offset_rad` を持たせています。

論文の Table 2 からは FSR と遅延長は分かりますが、簡易 Python モデルで必要な絶対位相基準までは完全には決まりません。実際の Lumerical モデルやレイアウトでは、MMI 内部位相、アームの共通長、曲げ・テーパー、ポート規約などにより自然に含まれる位相です。

この実装では、各 MZI 段について 0 から 2π までを探索し、Figure 1 のポート割当に沿って目標チャンネルが目的ポートに出るような位相オフセットを選びます。

つまり、

```text
dL               : 論文式・Table 2 に基づく周期設定
phase_offset_rad : 簡易モデルに不足している絶対位相基準の補正
```

という役割です。

## Output Files

`--csv` を指定すると、サンプリングしたスペクトルを CSV で保存します。

`--plot-file` を指定すると、matplotlib で PNG を保存します。

例:

```bash
uv run python main.py \
  --architecture lan-wdm \
  --span 45 \
  --points 451 \
  --csv out/lan_wdm_paper.csv \
  --plot-file out/lan_wdm_paper.png
```

## Notes

このプログラムは設計検討用の簡易モデルです。導波路分散の高次成分、MMI の波長依存、曲げ損失、モード変換、製造誤差、実際のレイアウト由来の位相は厳密には含んでいません。

実製造向けの最終検証には、FDTD/FDE/EME/INTERCONNECT などのフォトニック回路シミュレータでの確認が必要です。
