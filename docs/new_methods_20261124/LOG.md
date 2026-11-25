# LOG：ncRNA 新模型 / 新方法 deep research（2026-11-24 起）

- 主题短语（文件命名用）：`ncRNA-new-models`

- 作者 2026-11-24：「我并不在乎这个注释方法的路线怎么样，重点在于我们能不能提供一个新的模型、一种新的方法去研究 ncRNA。可以是序列、进化、测序数据、泛基因组、序列转图像……重点是找一种新的模型方法帮助学界更好地研究 ncRNA。」
- 流程：three-ai-deep-research skill。第 1 轮中立提示 `PROMPT_deep_research_EN.txt`（不含项目细节和已有结论），之后按第 2–5 轮讨论。
- 已尝试清单（第 4 轮才公开）：`09_deep_research_20261123/ATTEMPTS_SUMMARY_EN.md`，以及之后的 AUTH / E1 / IDENT / RNA_fate / 功能筛选 / 分层评价。

## 对话链接
- 第 1 轮 ChatGPT 深度研究（普通对话，非项目；作者 2026-11-24 决定保留，不在项目里重做）：https://chatgpt.com/c/6b060500-fe68-83eb-8308-ba8f55dde191
- 第 1 轮 Gemini Deep Research（模型「自动」，来源只勾搜索）：https://gemini.google.com/app/30e30f82d0db1a61
- 按作者 2026-11-24 新规则：不做 Claude 网页端调研、不做本地子代理调研（之前开的 4 个本地子代理在会话重启时中断，未产出，不补）。
- 第 1 轮 ChatGPT 报告：deep-research-report.md → deep_research/20261124_ncRNA-new-models_ChatGPT_R1.md（   53663 字节）
- 第 1 轮 Gemini 报告：页面 Blob 下载 ncRNA16_Gemini_R1_20261124-2342.md → deep_research/20261124_ncRNA-new-models_Gemini_R1.md（报告面板的 innerText；2340 版误取了提示面板，已删除，非原始 Markdown；表格为纯文本）
- 第 1 轮合并：`20261124_ncRNA-new-models_R1_merged.md`（关键事实 Europe PMC / GitHub API 核对；ChatGPT 导出的引用只有占位符）。
- 第 2 轮互评已发出（同一对话、普通模式）：提示 `discussion/PROMPT_round2_ChatGPT_EN.txt`、`PROMPT_round2_Gemini_EN.txt`；要求 ChatGPT 补 DOI/PMID，Gemini 换掉弱来源。
- 第 2 轮 ChatGPT 回复：Blob 下载被浏览器静默拦截，改为页面抽取，存要点存档 discussion/20261124_ncRNA-new-models_ChatGPT_R2.md（非逐字原文）。
- 第 2 轮 Gemini 回复：Deep Research 模式卡住 → 取消 + 刷新后显示；要点存档 discussion/20261124_ncRNA-new-models_Gemini_R2.md（非逐字）。主持人判断：discussion/20261124_ncRNA-new-models_Claude_R2.md。
- 第 3 轮（仓库 + 跨领域，候选 A–D）已发给两家：discussion/PROMPT_round3_EN.txt；给 Gemini 加了「不要显示内部起草过程」。
- 第 3 轮 Gemini：生成 20 多分钟后回复「身为一个语言模型，我没法提供这方面的帮助」（模型为自动，非 Deep Research 芯片）；同一对话精简重发一次。
- 第 3 轮 Gemini 重发后只输出 541 字符就中断。原 Deep Research 对话已不可靠（第 2 轮卡住、第 3 轮拒答再截断），判断是上下文过长。
- 改动（2026-11-25）：
  - ChatGPT 在原对话发第 4+5 轮合并（揭底 + 打分收敛），提示 `discussion/PROMPT_round45_ChatGPT_EN.txt`；
  - Gemini 开新普通对话发第 3–5 轮合并（附前情摘要 + 揭底），提示 `discussion/PROMPT_round345_Gemini_newchat_EN.txt`，链接 https://gemini.google.com/app/53b4b54c2f01f991。
- 揭底时核对：GitHub 仓库 JiWang-Unige/ncRNA 当前为 **public**（GitHub API `private: false`），而 ATTEMPTS_SUMMARY 里写的是 private。需告知作者，不擅自更改。
- 第 4+5 轮 ChatGPT（读了仓库）：B > A > C > D，要点存档 discussion/20261124_ncRNA-new-models_ChatGPT_R45.md。
- 第 3–5 轮 Gemini（新对话）：D/E > B > A，要点存档 discussion/20261125_ncRNA-new-models_Gemini_R345.md。
- 收束：20261125_ncRNA-new-models_SYNTHESIS.md。共识首选 B（区室反应性变化）；D 分歧，用 2 天审计判定；A 备选；C 搁置。待作者决定预研、下载和仓库可见性。
