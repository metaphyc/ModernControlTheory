# 现代控制理论作业解答

本仓库整理现代控制理论课程作业，参考本地《线性系统理论》教材，使用共享 LaTeX 模板按章排版。解答范围以 `assets/习题一.pdf` 至 `assets/习题七.pdf` 中红框标出的题号为准：框选整题时包含全部小问，框选小问时仅解答该小问。

## 作业索引

下表保留原题号。正文在作业标题之后直接开始第一题；选题及来源信息统一列在这里，PDF 页码从各题面文件的第 1 页开始计数。

| 章节 | 内容 | 红框选题 | 解答 |
| --- | --- | --- | --- |
| 习题一 | 状态空间描述及运动分析 | 2、4、7、10、15、16 | [PDF](homework/ch01/main.pdf) · [LaTeX](homework/ch01/main.tex) |
| 习题二 | 能控性 | 7（2）、9、11、14、15、20、22 | [PDF](homework/ch02/main.pdf) · [LaTeX](homework/ch02/main.tex) |
| 习题三 | 状态反馈与闭环极点配置 | 1、6、10、12、13 | [PDF](homework/ch03/main.pdf) · [LaTeX](homework/ch03/main.tex) |
| 习题四 | 能观性 | 6、9、10、13、17、18 | [PDF](homework/ch04/main.pdf) · [LaTeX](homework/ch04/main.tex) |
| 习题五 | 能控性、能观性与传递函数 | 1（c）、4、5、6（4）、6（7）、10 | [PDF](homework/ch05/main.pdf) · [LaTeX](homework/ch05/main.tex) |
| 习题六 | 状态观测器 | 1、2、6、7 | [PDF](homework/ch06/main.pdf) · [LaTeX](homework/ch06/main.tex) |
| 习题七 | 线性二次型最优控制 | 2、5 | [PDF](homework/ch07/main.pdf) · [LaTeX](homework/ch07/main.tex) |

习题五第 5 题包含（a）、（b）两种传递函数各自的四种实现；习题七第 2 题的性能指标续在题面第 2 页。解答均按这些完整范围整理。

| 题面文件 | 所选题目所在 PDF 页码 |
| --- | --- |
| 习题一.pdf | 第 1 页：2；第 2 页：4、7、10；第 3 页：15、16 |
| 习题二.pdf | 第 1 页：7（2）；第 2 页：9、11、14、15；第 3 页：20；第 4 页：22 |
| 习题三.pdf | 第 1 页：1；第 2 页：6；第 3 页：10、12、13 |
| 习题四.pdf | 第 1 页：6；第 2 页：9、10、13；第 3 页：17、18 |
| 习题五.pdf | 第 2 页：1（c）、4、5、6（4）；第 3 页：6（7）、10 |
| 习题六.pdf | 第 1 页：1、2；第 2 页：6、7 |
| 习题七.pdf | 第 1—2 页：2；第 2 页：5 |

## 目录与约定

```text
assets/                     原始题面 PDF（本地资料）
courseware/                 参考教材与课件（本地资料）
homework/
  _shared/preamble.tex      共享字体、页眉和题目/解答环境
  ch01/ ... ch07/
    main.tex                本章题意、推导及结果
    main.pdf                可直接阅读的解答
    verify.py               本章关键计算的符号核验
    figures/                本章题目所需的图形资源（如有）
PROMPT.md                   后续接手与新增作业的工作要求
```

原始 PDF 保留中文文件名。`assets/`、`courseware/` 的内容由 `.gitignore` 忽略，只跟踪 `.gitkeep` 以保留目录；解答 PDF 正常纳入版本管理，LaTeX 中间文件忽略。阅读或重新编译现有解答无需原始 PDF；核对题面或新增解答时需要本地资料。

采用与教材一致的正号状态反馈 `u = Kx + v`，闭环矩阵为 `A + BK`。坐标变换、观测器误差及最优控制的代价约定在各题中明确说明。解答按方法选择、分步推导、结论和核对组织，补全关键矩阵运算、系数匹配和证明依据，便于逐步阅读与复算。

题目中的图形应尽量忠实呈现，保留原布局、连接、箭头、符号及图号。习题一第 2 题的图 1.7 已从原题清晰裁取，资源位于 `homework/ch01/figures/fig-1-7.png`，编译时无需读取整份原始题面。

`homework/_shared/preamble.tex` 是用户原有模板，必须保持原样。后续作业禁止修改该文件，也禁止在章节或辅助文件中重定义、打补丁或覆盖其字体、页眉、题目框、解答框和标题样式。章节只使用模板已有接口填写题目与解答，保留原有固定文字和编号；不添加章节选题说明或统一约定的开场段落。

## 编译

需要 XeLaTeX 与 `latexmk`。TeX Live 完整安装包含模板用到的 `ctex`、Fandol 中文字体、`tcolorbox`、TikZ 和 `pgfplots` 等依赖。

编译单章：

```bash
cd homework/ch01
latexmk -xelatex -interaction=nonstopmode -halt-on-error main.tex
```

在仓库根目录编译全部七章：

```bash
for dir in homework/ch0[1-7]; do
  (cd "$dir" && latexmk -xelatex -interaction=nonstopmode -halt-on-error main.tex) || break
done
```

交付前检查编译错误、缺失字符、未定义引用和公式溢出，并查看生成 PDF 的排版。排版调整应限于正文的公式换行等内容，共享模板及其样式保持原样。

## 计算核验

各章 `verify.py` 使用 Python 3 与 SymPy，核对传递函数、状态响应、秩与变换、极点、观测器误差动态及 Riccati 方程等计算；证明题的论证保留在正文中。

在已安装 SymPy 的 Python 环境中，从仓库根目录执行：

```bash
for script in homework/ch0[1-7]/verify.py; do
  python3 "$script" || break
done
```

后续作业处理规则见 [PROMPT.md](PROMPT.md)。
