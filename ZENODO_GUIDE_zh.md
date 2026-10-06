# 你本人操作的 GitHub → Zenodo 归档

官方帮助核查日期：2026-10-05。GitHub 源码仓库：
<https://github.com/williamjay1/project40-reusable-transport-packaging>。
源码上传不等于正式 Release、Zenodo 归档或 DOI；本次不创建 Release、不操作 Zenodo。
用户于 2026-10-06 确认唯一的人类项目贡献者为 `junjiezhang`，保持此拼写；
不推断姓名拆分、ORCID、单位或 GitHub 账号。来源论文及数据提供者保留引用，
不列作本项目共同作者。原创材料公开许可尚待确认，当前未擅自授权 MIT/CC BY。

1. 检查上述仓库中的代码、数据及真实作者，确认原创材料许可，保留来源条款。
   用实际 commit 或正式 tag 记录所用版本；本地完整复现记录已随包提供。
   GitHub 仓库需满足 Zenodo 集成要求；组织仓库可能需组织管理员批准授权。
2. 你本人登录 Zenodo，连接 GitHub 账号；在个人菜单的 GitHub 页面`Sync now`，
   找到选定仓库并打开对应集成开关。
3. 检查 GitHub 默认分支的`CITATION.cff`，确保真实作者、题名、版本、仓库 URL
   和许可无误。现有生成工具要求已确认的姓名拆分；不能自行把 `junjiezhang`
   拆成 given/family names。当前确认记录不是可直接提交的 CFF 或 Zenodo 描述文件。
   此包不激活`.zenodo.json`，避免它覆盖 CFF；如你后来选择使用两者，
   Zenodo 官方说明：归档元数据只采用`.zenodo.json`，不会合并两份文件。
   作者元数据/许可修改后，运行`tools/refresh_manifest.py --confirm-metadata-only`
   及`tools/audit_repository.py`，保留此前 manifest；不要通过刷新清单隐藏模型改动。
4. 对已核验 commit 创建一个正式版本 tag/release。启用集成后，新 release 会由
   Zenodo 接收归档。不要把“提交代码”或“建立 tag”本身当成 DOI 已生成；确认实际
   归档记录/文件/许可/作者与该 release 一致。
5. 复制该版本实际获得的 DOI；同时保留 concept DOI（若平台提供）。论文复现声明
   应优先引用实际计算使用的具体版本 DOI，并列出 GitHub URL、tag/commit。不要用
   concept DOI 模糊哪个版本产生了正文结果。
6. 拿到 DOI 后，将真实值写回引用元数据与论文。若需改变已经归档的代码/数据，
   创建新版本并如实说明变更，不把新结果混入旧 DOI 所标识的版本。

获取 DOI 不证明数据可再分发、结论有效或期刊已接收。本仓库为 mixed-rights
软件/数据复现材料；原始文档和商业过程数据库未包含，第三方来源仍须引用。

官方入口：

- [启用仓库集成](https://help.zenodo.org/docs/github/enable-repository/)
- [Zenodo 软件元数据与优先规则](https://help.zenodo.org/docs/github/describe-software/)
- [GitHub 持久标识符/归档说明](https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content)
- [CITATION.cff 支持与校验](https://help.zenodo.org/docs/github/describe-software/citation-file/)

上述流程须由你本人进行账号连接、审批和归档确认；本包无账号 cookie、令牌、
上传密钥或自动发布命令。若集成界面改变，以当时官方帮助为准。
