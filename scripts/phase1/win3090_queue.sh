#!/bin/bash
# Start the v1.8 win3090 queue:  ssh win3090-wsl 'cd ~/ncrna && nohup bash scripts/phase1/win3090_queue.sh >> logs/formal_queue.log 2>&1 < /dev/null &'
cd ~/ncrna && source win3090.env && exec bash scripts/phase1/formal_queue.sh queues/win3090_v18.txt
