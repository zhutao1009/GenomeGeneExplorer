"""生成一个演示用的数据压缩包（``examples/demo_bundle.zip``）。

包里的文件名故意起得**毫无提示性**（``sequences_a.fna`` / ``sequences_b.fna``），
用来演示"导入压缩包后靠内容自动识别文件"这件事：程序会通过字母构成、序列名与
注释的对应关系、以及"记录长度是否等于注释里的 mRNA/CDS 长度"来判断每个文件。

用法：
    python tools/make_bundle.py            # 写到 examples/demo_bundle.zip
    python tools/make_bundle.py -o x.zip
"""

from __future__ import annotations

import argparse
import os
import sys
import zipfile

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

README = """演示数据集（Genome Gene Explorer）

本包用于试用「文件 → 导入压缩包（自动识别文件）」。

  assembly_v1.fna      基因组序列（FASTA，3 条序列，约 21 kb）
  annotation.gff3      基因注释（GFF3，11 个基因，含正负链 / 多转录本 / 边界情况）
  sequences_a.fna      参考 mRNA（故意没有写 mrna，靠内容识别）
  sequences_b.fna      参考 CDS（故意没有写 cds，靠长度与注释比对识别）
  peptides.faa         参考蛋白

导入后程序会列出每个文件被判定成什么、依据是什么，确认无误再点「开始导入」。
数据来自 examples/make_example_data.py 的确定性合成数据，可安全用于演示。
"""


def build(out_path: str) -> str:
    src = os.path.join(HERE, "examples")
    files = [
        ("assembly_v1.fna", "demo_genome.fa"),
        ("annotation.gff3", "demo_annotation.gff3"),
        ("sequences_a.fna", "demo_mrna.fa"),
        ("sequences_b.fna", "demo_cds.fa"),
        ("peptides.faa", "demo_protein.fa"),
    ]
    missing = [s for _n, s in files if not os.path.isfile(os.path.join(src, s))]
    if missing:
        raise SystemExit(
            "缺少示例数据：" + "、".join(missing) +
            "\n请先运行 python examples/make_example_data.py")

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with zipfile.ZipFile(out_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for member, source in files:
            zf.write(os.path.join(src, source), f"demo_dataset/{member}")
        zf.writestr("demo_dataset/README.txt", README)
        zf.writestr("__MACOSX/._junk", "macOS 打包残留，导入时会被忽略")
    return out_path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="生成演示用数据压缩包")
    ap.add_argument("-o", "--out", default=os.path.join(HERE, "examples", "demo_bundle.zip"))
    args = ap.parse_args(argv)
    path = build(args.out)
    size = os.path.getsize(path)
    print(f"已生成 {path}（{size / 1024:.1f} KB）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
