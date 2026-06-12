# vtkIsotropicRemeshingFilter

VTK 用のアイソトロピック（等方性）リメッシングフィルターです。

Botsch & Kobbelt, *"A Remeshing Approach to Multiresolution Modeling"*
(Symposium on Geometry Processing 2004) で発表されたインクリメンタル・
リメッシングアルゴリズムを、公開論文のみに基づいて VTK のデータ構造で
独立に実装したものです（外部ジオメトリ処理ライブラリへの依存・参照は
ありません。BSD ライセンス互換）。各反復で以下を行います。

1. **Split** — 目標エッジ長の 4/3 より長いエッジを中点で分割
2. **Collapse** — 目標エッジ長の 4/5 より短いエッジを中点へ縮約
   （4/3 を超えるエッジが生じる場合・多様体性が壊れる場合は抑止）
3. **Flip** — 頂点バレンスを 6（境界では 4）に近づけるエッジフリップ
4. **Tangential relaxation** — 一様ラプラシアンの接平面射影による平滑化
5. **Projection** — `vtkStaticCellLocator` による元サーフェスへの再投影

内部ではハーフエッジ構造を構築して局所操作（split / collapse / flip）を
行います。

## 特長・制約

- **境界の保護**: 境界エッジは常に拘束されます。境界頂点は境界ポリライン
  に沿ってのみ移動し、幾何学的なコーナー（折れ点）は固定されます。
- **フィーチャー保護**（オプション）: `ProtectFeaturesOn()` で、二面角が
  `FeatureAngle` を超えるエッジを境界と同様に保護します。
- 入力は **多様体かつ向きが一貫した** サーフェスメッシュである必要があり
  ます。非三角形セルは内部で三角形化され、重複点はマージされます。
  向きが不揃いな場合は事前に `vtkPolyDataNormals`（`ConsistencyOn`）を
  適用してください。
- **属性の補間**（`InterpolateAttributes`、既定ON）: 点データは出力頂点の
  元サーフェス上最近接点における重心座標補間で、セルデータは出力三角形
  重心に最も近い元セルからのコピーで引き継がれます（出力頂点は最終投影で
  元サーフェス上に乗るため、点データは実質的に正確な線形補間になります）。

## パラメーター

| パラメーター | 既定値 | 説明 |
|---|---|---|
| `TargetEdgeLength` | `0`（自動） | 目標エッジ長。`<= 0` の場合は入力の平均エッジ長を使用 |
| `NumberOfIterations` | `3` | リメッシング反復回数 |
| `NumberOfRelaxationSteps` | `1` | 反復ごとの接線平滑化回数 |
| `DoProject` | `true` | 元サーフェスへの再投影の有無 |
| `InterpolateAttributes` | `true` | 点データ・セルデータの出力への補間の有無 |
| `ProtectFeatures` | `false` | フィーチャーエッジ保護の有無 |
| `FeatureAngle` | `60`（度） | フィーチャー／コーナー判定の角度しきい値 |

## 使用例

```cpp
#include "vtkIsotropicRemeshingFilter.h"

vtkNew<vtkIsotropicRemeshingFilter> remesh;
remesh->SetInputData(surface);          // vtkPolyData
remesh->SetTargetEdgeLength(0.5);       // モデルのスケールに合わせて指定
remesh->SetNumberOfIterations(5);
remesh->ProtectFeaturesOn();            // CADモデル等のエッジを保持する場合
remesh->SetFeatureAngle(60.0);
remesh->Update();
vtkPolyData* result = remesh->GetOutput();
```

## ビルド

```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 -DVTK_DIR=<VTKのビルド/インストールディレクトリ>
cmake --build build --config Release
ctest --test-dir build -C Release
```

## コマンドラインツール

`RemeshSurface`（`BUILD_EXAMPLES=ON` でビルド）でファイルを直接リメッシュ
できます。

```
RemeshSurface input.{vtp|stl|ply|obj} output.{vtp|stl|ply} [targetEdgeLength] [iterations] [featureAngle]
```

- `targetEdgeLength` を省略（または 0）すると入力の平均エッジ長を使用
- `featureAngle` を指定するとフィーチャー保護が有効になります

## ディレクトリ構成

```
src/        vtkIsotropicRemeshingFilter (ライブラリ本体)
tests/      回帰テスト（球・境界付き平面・立方体のフィーチャー保護）
examples/   RemeshSurface コマンドラインツール
```

## 参考文献

- M. Botsch, L. Kobbelt, "A Remeshing Approach to Multiresolution Modeling",
  Symposium on Geometry Processing, 2004.
