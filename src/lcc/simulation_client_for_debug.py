# Copyright 2015 gRPC authors.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""The Python implementation of the GRPC helloworld.Greeter client."""

from __future__ import print_function

import asyncio
import logging

import sys
sys.path.insert(0, './protos')

import grpc
from grpc import aio
from protos import lccsimulation_pb2
from protos import lccsimulation_pb2_grpc

import time

ChannelName = 'localhost:50051'

async def sim_run() -> None:
    async with grpc.aio.insecure_channel(ChannelName) as channel:
        stub = lccsimulation_pb2_grpc.LCCSimulationStub(channel)
        response_iterator = stub.Run(lccsimulation_pb2.MEmpty())

        prevt = time.time()
        elapsed_time = 0
        count = 0
        async for response in response_iterator:
            curt = time.time()
            deltat = curt - prevt

            prevt = curt
            elapsed_time += deltat

            if elapsed_time >= 1:
                
                print("Received %d sim data in 1 second" % (count))

                elapsed_time = 0
                count = 0

            # for key, val in response.State.Objects.items():
            #     print(val.ID, val.Position)

            print('Simulation Step: ', response.State.SimulationStep, '\tSimulated Time', response.State.SimulatedTime)
            count += 1

def run():
    while(True):
        cmd = input("Server Function To Run: ")
        if cmd == "Initialize" or cmd == "0":
            Max_Message_Length = 10000000
            channel = grpc.insecure_channel(
                ChannelName,
                options=[
                    ('grpc.max_send_message_length', Max_Message_Length),
                    ('grpc.max_receive_message_length', Max_Message_Length),
                ],
            )
            stub = lccsimulation_pb2_grpc.LCCSimulationStub(channel)
            response = stub.Initialize(lccsimulation_pb2.MEmpty())
            print(response)

        elif cmd == "Run" or cmd == "1":
            asyncio.run(sim_run())

        elif cmd == "Pause" or cmd == "2":
            with grpc.insecure_channel(ChannelName) as channel:
                stub = lccsimulation_pb2_grpc.LCCSimulationStub(channel)
                response = stub.Pause(lccsimulation_pb2.MEmpty())

        elif cmd == "Stop" or cmd == "3":
            with grpc.insecure_channel(ChannelName) as channel:
                stub = lccsimulation_pb2_grpc.LCCSimulationStub(channel)
                response = stub.Stop(lccsimulation_pb2.MEmpty())

        elif cmd[:4] == "Plot":
            cmd_parsed = cmd.split(' ')
            Query = cmd_parsed[1]   # Currently, it does not accept indexing for the organism

            Max_Message_Length = 10000000
            channel = grpc.insecure_channel(
                ChannelName,
                options=[
                    ('grpc.max_send_message_length', Max_Message_Length),
                    ('grpc.max_receive_message_length', Max_Message_Length),
                ],
            )
            stub = lccsimulation_pb2_grpc.LCCSimulationStub(channel)
            response = stub.GetStaticPlotData(lccsimulation_pb2.MStaticPlotRequest(Identifier=Query))

        elif cmd[:5] == "Table":
            cmd_parsed = cmd.split(' ')
            Query = cmd_parsed[1]   # Currently, it does not accept indexing for the organism

            Max_Message_Length = 10000000
            channel = grpc.insecure_channel(
                ChannelName,
                options=[
                    ('grpc.max_send_message_length', Max_Message_Length),
                    ('grpc.max_receive_message_length', Max_Message_Length),
                ],
            )
            stub = lccsimulation_pb2_grpc.LCCSimulationStub(channel)
            response = stub.GetStaticTableData(lccsimulation_pb2.MStaticTableRequest(Identifier=Query))

        elif cmd == "Quit" or cmd == "-1":
            sys.exit(0)

        else:
            print("Function Not Found.")


if __name__ == '__main__':
    logging.basicConfig()
    run()
