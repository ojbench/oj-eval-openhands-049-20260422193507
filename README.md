# Problem 049 - 迅影的机器人（寒假 ver）

**ACMOJ Problem ID**: 2571

## Table of Contents

- [Problem 049 - 迅影的机器人（寒假 ver）](#problem-049-迅影的机器人寒假-ver)
  - [Table of Contents](#table-of-contents)
  - [Introduction](#introduction)
    - [Background](#background)
  - [Assignment Description](#assignment-description)
    - [Grade Composition](#grade-composition)
  - [Assignment Requirements](#assignment-requirements)
    - [Input Format](#input-format)
    - [Output Format](#output-format)
    - [Samples](#samples)
    - [Data Constraints](#data-constraints)
  - [Per-Testcase Resource Limits](#per-testcase-resource-limits)
  - [Test Data](#test-data)
  - [Submission Requirements](#submission-requirements)
    - [Submission Format](#submission-format)
    - [Submission Guidelines](#submission-guidelines)
    - [Evaluation Notes](#evaluation-notes)
    - [Academic Integrity](#academic-integrity)

## Introduction

### Background

Multi-robot collision avoidance (holiday version). Open-ended optimization problem.

## Assignment Description

### Grade Composition

| Grading Component | Percentage |
| :--: | :--: |
| Pass **2571. 迅影的机器人（寒假 ver）** | 80% |
| Code Review | 20% |

Here are several points that need clarification:

- In the Code Review, we will **strictly examine your code style and repository organization structure, etc.**

- This assignment provides some sample data for testing, stored in the `/workspace/data/049/data_test/` directory. Note that these are not the test cases on the Online Judge. Passing all local test cases does not guarantee that you will pass the OJ tests.

- Besides the provided sample data, we also encourage you to design your own test data based on your program logic to assist debugging.

## Assignment Requirements

### Problem Description

*本题为 PPCA 选题，你不应该期待满分。本题寒假版满分为 200 分。*

近日，迅影从黑洞余弦的庭院购入了一批小兔子机器人，用来帮助他处理政务。但在实际到手后他才发现，黑心的黑洞余弦并没有给这些机器人内置移动程序，它们现在无法移动。

为了避免国家的要务被耽搁，迅影现命令你为这一批小兔子机器人实现移动的逻辑，从而使得它们在不会相互碰撞的前提下能够各自到达想去的地方。

善良的黑洞余弦帮助你对问题进行了以下建模：
* 在二维平面内考虑，认为每个机器人都是圆形；
* 不考虑加速度，即速度是可以突变的；
* 每个机器人都具有内部属性和外部属性，**内部属性不可被其他机器人观测，外部属性则可被观测**；
* 内部属性包括机器人自身可以达到的最大速度 $v^{max}$ 和机器人的目标位置坐标 $\vec{r}^{tar}$；
* 外部属性为机器人的当前位置坐标 $\vec{r}^{cur}$ 、机器人的半径  $R$ 、机器人的当前速度 $\vec{v}^{cur}$，以及上个时刻机器人的碰撞/超速信息；
* 总机器人数量是全局可知的；
* **目标：每隔时间间隔 $\tau$ 后改变机器人的速度，使得接下来的 $\tau$ 之内，每个机器人都不与其他机器人相撞、且速度不会超过该机器人的最大速度限制，同时在最大轮数限制内使得所有的机器人均到达目标位置。**

### 任务目标与说明

体贴的黑洞余弦向你提供了一个框架(见[attachment]第三次机考B卷下发文件.zip[/attachment] / [jbox](https://jbox.sjtu.edu.cn/l/y1gn6z "下发文件")中的`src`文件夹)，你需要实现的是
`src.hpp` 中的 `Controller::get_v_next()` 函数，即根据当前机器人的状态和其他机器人的外部属性状态，**在每个时刻决定接下来的时间间隔 $\tau$ 内该机器人的速度**。

具体来说，在每个时刻开始前，程序会按照编号**从小到大**的顺序，对每个机器人调用 `Controller::get_v_next()` 函数以获得接下来的 $\tau$ 内该机器人的速度。然后程序会进行超速检测和碰撞检测，**如果在接下来的 $\tau$ 内所有机器人均未发生超速和碰撞，那么程序会更新各个机器人的状态；否则，该步无效，所有机器人的状态均不会改变**。机器人在到达目标地点后并不会主动停下来，而是会继续模拟直至所有机器人均到达或者达到最大时间限制。

框架向你提供了一个向量类 `Vec` (在`math.h`中)，实现了二维向量的加减、数乘、点乘、叉乘、求模、旋转等操作。你需要了解 `Vec` 类并使用它。

在 `Controller` 类中，你可以通过**访问它的成员变量来获取当前机器人的属性状态**(上个时刻结束时的位置 `pos_cur`，上个时刻的位置 `v_cur`，目标位置 `pos_tar`，半径 `r`，最大速度 `v_max`，以及机器人编号 `id` (0-based))。对于其他机器人的外部属性状态，你可以通过类成员中的 `*monitor` 来访问，这是一个全局的沟通渠道。黑洞余弦向你提供了以下访问接口：

```cpp
bool Monitor::get_speeding(int id) const; // 若机器人id超速，返回true
std::vector<int> Monitor::get_collision(int id) const; // 返回在上个时刻的模拟中和机器人id相碰的机器人的编号(若无，返回空vector)
bool Monitor::get_warning() const; // 若上个时刻机器人id的速度合法，返回true
Vec Monitor::get_pos_cur(int id) const; // 返回机器人id在上个时刻结束时的位置
Vec Monitor::get_v_cur(int id) const; // 返回机器人id在上个时刻的速度
double Monitor::get_r(int id) const; // 返回机器人id的半径
bool Monitor::get_done() const; // 如果全部机器人均已达到目标位置，返回true
int Monitor::get_robot_number() const; // 返回机器人总数
int Monitor::get_test_id() const; // 返回当前的测试点编号
```

为了防止你走投无路，可爱的黑洞余弦还额外向你提供了一个**64bit的全局数据共享区**。你可以调用 `Monitor::write_share_info()` 和 `Monitor::use_share_info()` 来分别对共享区域的数据进行读写。**但是注意，每次使用共享区都会进行一定的扣分，具体细则将在下文说明。**

另外，为了防止你一头雾水，在迅影的劝说下，黑洞余弦还在[下发文件](https://jbox.sjtu.edu.cn/l/y1gn6z "下发文件")中给了你**部分测试点**。你可以在下发文件中的`testcases`文件夹中找到它们。希望这些下发测试点会对你有用！

【**Tips1：除了 `src.hpp` 中的指定部分，禁止修改其他的任何代码！禁止使用任何的静态变量或全局变量来试图实现不同机器人之间的数据共享！禁止输出任何内容来尝试攻击测评程序！所有的非法行为均会被助教人工检查并void。**】

【Tips2：请勿忘记时间间隔 $\tau$，你计算出的速度应该是 $ds/dt$。在本题中，$\tau$ 在 `config.h` 中，被规定为 `TIME_INTERVAL = 0.1`。】

【Tips3：请务必关注评分标准公式并仔细研究，以期实现最优得分策略。】

【Tips4：请务必充分利用下发数据。你可以通过 `Monitor` 类知道当前测试点编号。】

【Tips5：关于浮点数，评测程序均已经加入一定的、合理的误差容忍范围。另外你可能需要特别注意 `inf`、`nan` 等特殊情况。】


### 可视化调试

为了本地调试方便，我们另外特地提供了一个基于python的本地可视化途径。

在程序成功运行之后，会在当前目录的文件夹中生成`trace_info`，它是你的机器人的运动数据。你可以直接根据它来进行调试(具体文件含义见`visualization/file.py`)，但更建议你运行`visualization/main.py`，从而生成一张`trace.gif`进行辅助调试。如果你觉得图片生成的太慢，可以更改`config.py`中的参数从而调整生成速度。【如果你使用CMake，`trace_info`可能会生成在`cmake-build-debug`文件夹中，这时你可能需要在`config.h` 或者 `visualization/config.py` 中更改文件路径】

如果你已经配置过python环境，直接按照你喜欢的方式运行`main.py`即可；如果你还没有配置过python环境，请在wsl中输入下列指令：
```shell
sudo apt install python3-pip
pip insatll matplotlib
```
完成后，输入
```bash
cd visualization
python3 main.py
```

即可运行。

**环境配置并不是我们所想要考察的内容，如果在环境配置上有任何问题，请及时与助教联系。**

### 评分标准

对于每个测试点，若其满分为 $S$，共有 $N$ 个机器人，最大轮数限制为 $T$，程序在 $t$ 轮后结束，结束时有 $n$ 个机器人到达了目标，共享区使用(读/写)了 $m$ 次。那么最终得分 $x$ 为:

$x=S \cdot \max\left\\{\left(0.7 \cdot \left(\frac{n}{N}\right)^2 + 0.3\cdot\left(1-\sqrt{\frac{t}{T}}\right)- \frac{m}{100} \cdot\left(0.05m + 4.95\right)\right),0\right\\}$

### Input Format

请按上述要求完成下发文件中的 src.hpp 并提交。

另外对下发样例的输入格式进行说明：

- 第一行一个正整数 $id$，表示测试点编号；

- 第二行一个正整数 $T$，表示最大运行轮数；

- 第三行一个正整数 $N$，表示机器人数量；

- 接下来的 $N$ 行，每行六个**浮点数** $x_{init}^i, y_{init}^i, x_{tar}^i,y_{tar}^i,r,v_{max}$，表示每个机器人的初始坐标，目标坐标，半径，最大速度。

### Output Format

|  测试组 | 测试点  | 分数 | 描述 |
| :------------: | :------------: | :--:| :--:|
|  简单 |  1-8 |  每个点 3.58 分| 最多只有两个机器人 |
|  普通 | 9-16  |  每个点 7.14 分 |最多只有六个机器人|
| 困难 | 17-24 | 每个点 14.29 分 | 不会超过十二个机器人|

保证所有数据点的输入均合法，保证所有数据点对应的机器人均可以在符合要求的前提下到达目标位置。

对于所有数据点，保证 $-1000 \leq x_{init}^i, y_{init}^i, x_{tar}^i,y_{tar}^i,r,v_{max} \leq 1000$。

### Samples

No sample data provided for this problem.

### Data Constraints

See the problem description for constraints.

## Per-Testcase Resource Limits

- **Time Limit (per test case)**: 15000 ms
- **Memory Limit (per test case)**: 256 MiB
- **Disk Usage**: No disk usage is permitted.

## Test Data

The test data for this problem is located at `/workspace/data/049/data_test/`.

Attachments and starter files for this problem are located at `/workspace/data/049/`.

## Submission Requirements

### Submission Format

Submit a single C++ source file to the Online Judge. The OJ provides its own `main.cpp` which `#include "src.hpp"`. Your submission will be written to `src.hpp` and compiled together with the OJ's driver code.

### Submission Guidelines

- The submitted code must be able to compile successfully through the above compilation process
- The program needs to be able to read data from standard input and write results to standard output
- Please ensure the code runs correctly within the given time and space limits
- **You must use C++ or C language** to implement this assignment

### Evaluation Notes

- The evaluation system will test your program using the provided test data
- The program output must exactly match the expected output (including format)
- Exceeding time or memory limits will be judged as the corresponding error type
- Please pay attention to the overall time performance of your code and the time complexity of each part of your algorithm.

### Academic Integrity

If any violations are found during evaluation or code review (including but not limited to using unconventional methods to pass test cases), your final score may be significantly reduced or become **0 points**.
