# subevent_motion_alg 主机构建 + 回归测试（与 Python 参考实现逐点对比）
#   make test        # 用 Python 导出测试向量，分别以 PAIR_ORDER = 1 / 2 编译并运行
# test/port/ 下是工程头文件（common_util.h、securec.h）的主机测试桩，集成到工程时使用工程自己的头文件。
CC      ?= gcc
CFLAGS  ?= -std=c99 -O2 -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wno-sign-conversion
PYTHON  ?= python3
SRCS    := subevent_motion_alg.c test/subevent_motion_test.c
DEPS    := $(SRCS) subevent_motion_alg.h test/port/common_util.h test/port/securec.h

all: subevent_motion_test_o1 subevent_motion_test_o2 replay_dump

# 用固件 [jdps_dump] 输入回放：./replay_dump firmware.log
replay_dump: subevent_motion_alg.c subevent_motion_alg.h test/replay_dump_main.c test/port/common_util.h test/port/securec.h
	$(CC) $(CFLAGS) -D_DEFAULT_SOURCE -I. -Itest/port -o $@ subevent_motion_alg.c test/replay_dump_main.c -lm

subevent_motion_test_o%: $(DEPS)
	$(CC) $(CFLAGS) -DSUBEVENT_MOTION_PAIR_ORDER=$* -D_DEFAULT_SOURCE -I. -Itest/port -o $@ $(SRCS) -lm

test/vectors.txt: ../cs_agc/algorithms.py ../cs_agc/export_vectors.py
	cd .. && $(PYTHON) -m cs_agc.export_vectors c/test/vectors.txt

test: all test/vectors.txt
	./subevent_motion_test_o1 test/vectors.txt
	./subevent_motion_test_o2 test/vectors.txt

clean:
	rm -f subevent_motion_test_o1 subevent_motion_test_o2 replay_dump test/vectors.txt

.PHONY: all test clean
